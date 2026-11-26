"""Vacancy-mediated dopant tracer diffusion on a periodic diamond lattice.

A small 3D diamond-cubic network is generated explicitly (four nearest
neighbours per site). One vacancy and one dopant tracer are evolved with
residence-time KMC. This is a pedagogical single-vacancy model; quantitative
Si:B / Si:P predictions require validated defect species, charge states,
Fermi-level effects, binding energies, and calibrated transition barriers.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result

K_B_EV = 8.617333262145e-5
_FCC = ((0, 0, 0), (0, 2, 2), (2, 0, 2), (2, 2, 0))
_A_TO_B = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))


def build_diamond_lattice(n_cells: int, lattice_constant_nm: float):
    if n_cells < 2:
        raise ValueError("n_cells must be >= 2")
    L = 4 * n_cells
    coords = set()
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                cell = (4 * i, 4 * j, 4 * k)
                for f in _FCC:
                    a = tuple((cell[d] + f[d]) % L for d in range(3))
                    b = tuple((a[d] + 1) % L for d in range(3))
                    coords.add(a)
                    coords.add(b)
    ordered = sorted(coords)
    index = {c: n for n, c in enumerate(ordered)}
    neighbours = np.empty((len(ordered), 4), dtype=np.int32)
    deltas = np.empty((len(ordered), 4, 3), dtype=np.int32)
    for idx, c in enumerate(ordered):
        is_a = (c[0] % 2 == 0)  # A sublattice coordinates are all even.
        steps = _A_TO_B if is_a else tuple(tuple(-x for x in d) for d in _A_TO_B)
        for m, d in enumerate(steps):
            cc = tuple((c[q] + d[q]) % L for q in range(3))
            if cc not in index:
                raise RuntimeError(f"diamond lattice neighbour missing: {cc}")
            neighbours[idx, m] = index[cc]
            deltas[idx, m] = d
    if np.any(np.diff(np.sort(neighbours, axis=1), axis=1) == 0):
        raise RuntimeError("duplicate nearest neighbours detected")
    # Every site must have four neighbours and the graph must be reciprocal.
    for i in range(len(ordered)):
        for j in neighbours[i]:
            if i not in neighbours[j]:
                raise RuntimeError("diamond neighbour graph is not reciprocal")
    pos_nm = np.asarray(ordered, dtype=float) * (lattice_constant_nm / 4.0)
    return np.asarray(ordered, dtype=np.int32), neighbours, deltas, pos_nm, L


def _run_one(neighbours, deltas, coords, L, a_nm, temp, em_v, em_d, nu, n_events, rng):
    n_sites = len(coords)
    vacancy = int(rng.integers(n_sites))
    dopant = int(rng.integers(n_sites - 1))
    if dopant >= vacancy:
        dopant += 1
    dopant_unwrapped = np.asarray(coords[dopant], dtype=float)
    t = 0.0
    sample_every = max(1, n_events // 120)
    times, msd, path = [0.0], [0.0], [dopant_unwrapped.copy()]
    rate_v = nu * np.exp(-em_v / (K_B_EV * temp))
    rate_d = nu * np.exp(-em_d / (K_B_EV * temp))
    for step in range(1, n_events + 1):
        dests = neighbours[vacancy]
        rates = np.where(dests == dopant, rate_d, rate_v)
        total = float(np.sum(rates))
        if total <= 0 or not np.isfinite(total):
            raise RuntimeError("non-positive or invalid total KMC rate")
        u = rng.random() * total
        choice = int(np.searchsorted(np.cumsum(rates), u, side="right"))
        choice = min(choice, 3)
        new_v = int(dests[choice])
        delta = deltas[vacancy, choice].astype(float)
        t += -np.log(max(rng.random(), np.finfo(float).tiny)) / total
        old_v = vacancy
        if new_v == dopant:
            # The dopant exchanges with the vacancy; its physical step is -delta.
            dopant_unwrapped -= delta
            dopant = old_v
        vacancy = new_v
        if step % sample_every == 0 or step == n_events:
            # Reference position is stored in path[0].
            dr = (dopant_unwrapped - path[0]) * (a_nm / 4.0)
            times.append(t)
            msd.append(float(np.dot(dr, dr)))
            path.append(dopant_unwrapped.copy())
    return np.asarray(times), np.asarray(msd), np.asarray(path), t


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2025)))
    n_cells = int(p.get("n_cells", 3))
    a_nm = float(p.get("lattice_constant_nm", 0.5431))
    temp = float(p.get("temperature_K", 1100.0))
    em_v = float(p.get("vacancy_migration_eV", 0.50))
    em_d = float(p.get("dopant_exchange_eV", 0.55))
    nu = float(p.get("attempt_frequency_s", 1e13))
    n_events = int(p.get("events_per_run", 5000))
    reps = int(p.get("replicates", 8))
    if temp <= 0 or a_nm <= 0 or em_v < 0 or em_d < 0 or nu <= 0 or n_events < 100 or reps < 2:
        raise ValueError("Invalid KMC parameters; use T>0, nonnegative barriers, >=100 events and >=2 replicates")
    coords, neigh, deltas, pos_nm, L = build_diamond_lattice(n_cells, a_nm)

    runs = []
    duration = []
    diffusion = []
    representative_path = None
    for rep in range(reps):
        trng = np.random.default_rng(int(rng.integers(0, 2**32 - 1)))
        times, msd, path, total_t = _run_one(neigh, deltas, coords, L, a_nm, temp, em_v, em_d, nu, n_events, trng)
        # Interpolate each replicate to a common normalized-time grid for ensemble averaging.
        tau = np.linspace(0.0, 1.0, 100)
        scaled_t = times / max(times[-1], np.finfo(float).tiny)
        msd_interp = np.interp(tau, scaled_t, msd)
        runs.append((tau, msd_interp, times, msd, path))
        duration.append(total_t)
        # Fit the late-time part of the tracer MSD; D = slope / 6 in 3D.
        start = max(2, int(0.25 * len(times)))
        slope = float(np.polyfit(times[start:], msd[start:], 1)[0]) if len(times[start:]) >= 2 else 0.0
        diffusion.append(slope / 6.0)
        if representative_path is None:
            representative_path = path

    tau = runs[0][0]
    mean_msd = np.mean([r[1] for r in runs], axis=0)
    mean_duration = float(np.mean(duration))
    mean_D = float(np.mean(diffusion))
    std_D = float(np.std(diffusion, ddof=1))
    mean_times = tau * mean_duration
    # Add a representative unwrapped tracer path for visualization and reproducibility.
    rep_times = runs[0][2]
    rep_path = runs[0][4]
    rep_tau = rep_times / max(rep_times[-1], np.finfo(float).tiny)
    rep_disp_nm = (rep_path - rep_path[0]) * (a_nm / 4.0)
    path_x = np.interp(tau, rep_tau, rep_disp_nm[:, 0])
    path_y = np.interp(tau, rep_tau, rep_disp_nm[:, 1])
    path_z = np.interp(tau, rep_tau, rep_disp_nm[:, 2])
    df = pd.DataFrame({"normalized_time": tau, "mean_physical_time_s": mean_times,
                       "ensemble_mean_tracer_MSD_nm2": mean_msd,
                       "representative_tracer_dx_nm": path_x,
                       "representative_tracer_dy_nm": path_y,
                       "representative_tracer_dz_nm": path_z})

    before = np.zeros((L, L), dtype=float)
    after = np.zeros((L, L), dtype=float)
    # Project all replicate tracer paths onto x-y to make a useful before/after map.
    for run in runs:
        path_i = run[4]
        start_xyz = path_i[0]
        before[int(start_xyz[0]) % L, int(start_xyz[1]) % L] += 1.0
        for point in path_i:
            after[int(point[0]) % L, int(point[1]) % L] += 1.0
    series = {
        "Ensemble tracer MSD": {"x": mean_times, "y": mean_msd,
                                "xlabel": "Physical time (s)", "ylabel": "Mean-square displacement (nm²)"},
        "Representative tracer x displacement": {"x": mean_times, "y": path_x,
                                "xlabel": "Physical time (s)", "ylabel": "Tracer displacement x (nm)"},
        "Representative tracer z displacement": {"x": mean_times, "y": path_z,
                                "xlabel": "Physical time (s)", "ylabel": "Tracer displacement z (nm)"},
    }
    metrics = {
        "mean_tracer_diffusivity_nm2_s": mean_D,
        "mean_tracer_diffusivity_m2_s": mean_D * 1e-18,
        "replicate_standard_deviation_nm2_s": std_D,
        "mean_simulated_duration_s": mean_duration,
        "diamond_lattice_sites": int(len(coords)),
        "coordination_number": 4,
        "replicates": reps,
        "events_per_run": n_events,
    }
    return make_result(
        name="dopant_diffusion", title="Vacancy-mediated dopant tracer diffusion in silicon",
        method="Kinetic Monte Carlo · diamond-cubic nearest-neighbour network",
        description="Tracks a dopant tracer exchanging with a mobile vacancy on a periodic diamond-cubic network.",
        before=before, after=after, data=df, series=series, metrics=metrics,
        visual_kind="lattice", before_label="Initial tracer position (2D projection)",
        after_label="Final tracer position (2D projection)",
        notes=[
            "The topology is diamond-cubic with four nearest neighbours per site; the model contains one vacancy and one dopant tracer per replicate.",
            "The reported diffusivity is a stochastic finite-run tracer estimate from the Einstein relation; uncertainty is estimated across independent runs.",
            "A single vacancy is not a realistic equilibrium vacancy concentration. Do not interpret this as a calibrated boron or phosphorus diffusivity without a defect-chemistry model.",
            "The model omits charge states, Fermi-level dependence, dopant-vacancy binding complexes, interstitial-mediated diffusion and concentration-dependent effects.",
        ],
        warnings=["Migration and exchange barriers are user-supplied placeholders unless sourced/calibrated for the specific dopant, defect charge state and temperature."],
        metadata={"temperature_K": temp, "vacancy_migration_eV": em_v,
                  "dopant_exchange_eV": em_d, "lattice_constant_nm": a_nm,
                  "attempt_frequency_s": nu, "seed": int(p.get("seed", 2025))},
    )
