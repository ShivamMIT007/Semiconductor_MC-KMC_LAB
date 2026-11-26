"""Coarse-grained object KMC for vacancy-interstitial recombination in a 2D lattice."""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result, periodic_neighbors_2d
K_B_EV = 8.617333262145e-5


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2029)))
    size = int(p.get("size", 36))
    pairs = int(p.get("initial_pairs", 45))
    n_events = int(p.get("events", 5000))
    temp = float(p.get("temperature_K", 900.0))
    em_v = float(p.get("vacancy_migration_eV", 0.65))
    em_i = float(p.get("interstitial_migration_eV", 0.25))
    er = float(p.get("recombination_barrier_eV", 0.05))
    nu = float(p.get("attempt_frequency_s", 1e12))
    nu_rec = float(p.get("recombination_prefactor_s", 1e12))
    if size < 6 or pairs < 1 or 2 * pairs >= size * size or n_events < 100 or temp <= 0 or min(em_v, em_i, er) < 0 or min(nu, nu_rec) <= 0:
        raise ValueError("Invalid defect KMC parameters")
    # 0 empty, 1 vacancy object, 2 self-interstitial object.
    state = np.zeros((size, size), dtype=np.int8)
    flat_sites = rng.choice(size * size, size=2 * pairs, replace=False)
    state.flat[flat_sites[:pairs]] = 1
    state.flat[flat_sites[pairs:]] = 2
    before = state.copy()
    rate_v = nu * np.exp(-em_v / (K_B_EV * temp))
    rate_i = nu * np.exp(-em_i / (K_B_EV * temp))
    rate_rec = nu_rec * np.exp(-er / (K_B_EV * temp))
    t = 0.0
    recombinations = 0
    sample_every = max(1, n_events // 100)
    records = []
    snapshots = [state.copy()]
    executed = 0
    for step in range(1, n_events + 1):
        events, rates = [], []
        for i, j in np.argwhere(state != 0):
            species = int(state[i, j])
            hop_rate = rate_v if species == 1 else rate_i
            for ni, nj in periodic_neighbors_2d(int(i), int(j), size, size):
                target = int(state[ni, nj])
                if target == 0:
                    events.append(("hop", int(i), int(j), ni, nj, species))
                    rates.append(hop_rate)
                elif species == 1 and target == 2:
                    events.append(("recombine", int(i), int(j), ni, nj, species))
                    rates.append(rate_rec)
        if not events:
            break
        total = float(np.sum(rates))
        idx = int(np.searchsorted(np.cumsum(rates), rng.random() * total, side="right"))
        idx = min(idx, len(events) - 1)
        t += -np.log(max(rng.random(), np.finfo(float).tiny)) / total
        kind, i, j, ni, nj, species = events[idx]
        if kind == "hop":
            state[ni, nj] = species
            state[i, j] = 0
        else:
            state[i, j] = 0
            state[ni, nj] = 0
            recombinations += 1
        executed += 1
        if step % sample_every == 0 or step == n_events:
            records.append({"event": step, "time_s": t,
                            "vacancy_count": int(np.count_nonzero(state == 1)),
                            "interstitial_count": int(np.count_nonzero(state == 2)),
                            "surviving_defects": int(np.count_nonzero(state)),
                            "recombination_events": recombinations})
            if step in {max(1, n_events // 4), max(1, n_events // 2), max(1, 3 * n_events // 4)}:
                snapshots.append(state.copy())
    df = pd.DataFrame(records)
    if not df.empty:
        series = {
            "Surviving defects": {"x": df["time_s"].to_numpy(), "y": df["surviving_defects"].to_numpy(),
                                  "xlabel": "Physical time (s)", "ylabel": "Defect objects"},
            "Vacancies and interstitials": {"x": df["time_s"].to_numpy(), "y": df["vacancy_count"].to_numpy(),
                                            "xlabel": "Physical time (s)", "ylabel": "Vacancy count"},
        }
    else:
        series = {}
    snapshots.append(state.copy())
    metrics = {"initial_defect_objects": int(np.count_nonzero(before)),
               "final_vacancies": int(np.count_nonzero(state == 1)),
               "final_interstitials": int(np.count_nonzero(state == 2)),
               "final_surviving_defects": int(np.count_nonzero(state)),
               "recombination_events": recombinations,
               "executed_KMC_events": executed,
               "physical_time_s": float(t)}
    return make_result(
        name="defect_annealing", title="Vacancy–interstitial defect annealing",
        method="Object KMC · defect migration and recombination",
        description="Tracks vacancy and self-interstitial objects on a periodic 2D coarse lattice.",
        before=before, after=state, data=df, series=series, metrics=metrics,
        visual_kind="defects", before_label="Initial defect population", after_label="After annealing events",
        snapshots=snapshots,
        notes=[
            "Vacancies and interstitials migrate with Arrhenius event rates; neighbouring V-I pairs can recombine.",
            "This is a 2D coarse-grained reaction-diffusion KMC demonstrator, not a silicon defect-cluster model.",
            "A research model must include species-specific migration, binding, charge states, cluster sizes, dissociation and appropriate 3D lattice geometry.",
        ],
        warnings=["Default barriers and prefactors are placeholders. Calibrate them using atomistic calculations and defect-annealing measurements."],
        metadata={"size": size, "initial_pairs": pairs, "temperature_K": temp,
                  "vacancy_migration_eV": em_v, "interstitial_migration_eV": em_i,
                  "recombination_barrier_eV": er, "seed": int(p.get("seed", 2029))},
    )
