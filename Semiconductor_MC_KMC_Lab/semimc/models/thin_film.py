"""Kinetic Monte Carlo solid-on-solid (SOS) thin-film deposition model.

A coarse-grained periodic square lattice evolves through deposition, surface
hopping and optional desorption events. It is useful for teaching rate-based
surface growth; it is not a chemically calibrated PVD/CVD/ALD process model.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result
from semimc.common import periodic_neighbors_2d
K_B_EV = 8.617333262145e-5


def _roughness(h):
    return float(np.std(h.astype(float)))


def _build_events(h, flux, nu, temp, e_diff, e_bond, e_des):
    nx, ny = h.shape
    events = []
    rates = []
    # Deposition arrival rate is flux ML/s per column.
    for i in range(nx):
        for j in range(ny):
            events.append(("deposit", i, j, -1, -1))
            rates.append(flux)
    # Surface diffusion: move one atom from the top of source column to neighbour.
    for i in range(nx):
        for j in range(ny):
            if h[i, j] <= 0:
                continue
            nbrs = periodic_neighbors_2d(i, j, nx, ny)
            lateral_bonds = sum(1 for ni, nj in nbrs if h[ni, nj] >= h[i, j])
            hop_rate = nu * np.exp(-(e_diff + lateral_bonds * e_bond) / (K_B_EV * temp))
            des_rate = nu * np.exp(-e_des / (K_B_EV * temp))
            for ni, nj in nbrs:
                events.append(("hop", i, j, ni, nj))
                rates.append(hop_rate)
            events.append(("desorb", i, j, -1, -1))
            rates.append(des_rate)
    return events, np.asarray(rates, dtype=float)


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2026)))
    size = int(p.get("size", 18))
    n_events = int(p.get("events", 6000))
    temp = float(p.get("temperature_K", 650.0))
    flux = float(p.get("flux_ml_s", 0.25))
    nu = float(p.get("effective_attempt_frequency_s", 1e5))
    e_diff = float(p.get("diffusion_barrier_eV", 0.65))
    e_bond = float(p.get("bond_penalty_eV", 0.08))
    e_des = float(p.get("desorption_barrier_eV", 1.8))
    if size < 4 or n_events < 100 or temp <= 0 or flux < 0 or nu <= 0 or min(e_diff, e_bond, e_des) < 0:
        raise ValueError("Invalid deposition parameters")
    h = np.zeros((size, size), dtype=np.int32)
    initial = h.copy()
    total_time = 0.0
    records = []
    snapshots = [h.copy()]
    sample_every = max(1, n_events // 100)
    deposition_count = 0
    hop_count = 0
    des_count = 0
    for step in range(1, n_events + 1):
        events, rates = _build_events(h, flux, nu, temp, e_diff, e_bond, e_des)
        total_rate = float(rates.sum())
        if not np.isfinite(total_rate) or total_rate <= 0:
            break
        r = rng.random() * total_rate
        idx = int(np.searchsorted(np.cumsum(rates), r, side="right"))
        idx = min(idx, len(events) - 1)
        event = events[idx]
        total_time += -np.log(max(rng.random(), np.finfo(float).tiny)) / total_rate
        kind, i, j, ni, nj = event
        if kind == "deposit":
            h[i, j] += 1
            deposition_count += 1
        elif kind == "hop":
            if h[i, j] > 0:
                h[i, j] -= 1
                h[ni, nj] += 1
                hop_count += 1
        elif kind == "desorb":
            if h[i, j] > 0:
                h[i, j] -= 1
                des_count += 1
        if step % sample_every == 0 or step == n_events:
            records.append({"event": step, "time_s": total_time,
                            "mean_thickness_ML": float(h.mean()),
                            "surface_roughness_ML": _roughness(h),
                            "max_height_ML": int(h.max()),
                            "deposition_events": deposition_count,
                            "surface_hops": hop_count,
                            "desorption_events": des_count})
            snapshots.append(h.copy())
    df = pd.DataFrame(records)
    series = {
        "Mean film thickness": {"x": df["time_s"].to_numpy(), "y": df["mean_thickness_ML"].to_numpy(),
                                "xlabel": "Physical time (s)", "ylabel": "Mean thickness (ML)"},
        "Surface roughness": {"x": df["time_s"].to_numpy(), "y": df["surface_roughness_ML"].to_numpy(),
                              "xlabel": "Physical time (s)", "ylabel": "RMS roughness (ML)"},
    }
    metrics = {
        "final_mean_thickness_ML": float(h.mean()),
        "final_surface_roughness_ML": _roughness(h),
        "final_max_height_ML": int(h.max()),
        "physical_simulation_time_s": float(total_time),
        "deposition_events": deposition_count,
        "surface_hops": hop_count,
        "desorption_events": des_count,
        "net_atom_count_in_columns": int(h.sum()),
    }
    return make_result(
        name="thin_film_deposition", title="KMC thin-film deposition (solid-on-solid)",
        method="Kinetic Monte Carlo · deposition / surface hopping / desorption",
        description="Coarse-grained surface growth on a periodic 2D solid-on-solid lattice.",
        before=initial, after=h, data=df, series=series, metrics=metrics,
        visual_kind="height", before_label="Bare substrate", after_label="Deposited film surface",
        snapshots=snapshots,
        notes=[
            "Deposition, hopping and desorption are explicit stochastic events with continuous KMC time advancement.",
            "The attempt frequency is an effective coarse-grained rate in the default example. Replace it with physically justified prefactors before quantitative comparison.",
            "This SOS model prohibits overhangs and does not explicitly represent precursor chemistry, gas transport, plasma species or crystallographic reconstruction.",
        ],
        warnings=["Use this as a surface-growth algorithm prototype. Real PVD/CVD/ALD requires process-specific fluxes, sticking probabilities, reaction networks and calibrated barriers."],
        metadata={"temperature_K": temp, "flux_ml_s": flux, "effective_attempt_frequency_s": nu,
                  "diffusion_barrier_eV": e_diff, "bond_penalty_eV": e_bond,
                  "desorption_barrier_eV": e_des, "seed": int(p.get("seed", 2026))},
    )
