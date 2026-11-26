"""2D Potts-model Monte Carlo grain growth (Metropolis acceptance).

One Monte Carlo step (MCS) means N attempted site updates. MCS is an algorithmic
clock, not automatically seconds. A calibrated mobility law is needed to map
MCS to physical time.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result, boundary_fraction, periodic_neighbors_2d


def _voronoi_initial(size: int, n_grains: int, rng: np.random.Generator) -> np.ndarray:
    seeds = rng.integers(0, size, size=(n_grains, 2))
    yy, xx = np.indices((size, size))
    dist = np.full((size, size), np.inf)
    lab = np.zeros((size, size), dtype=np.int32)
    for k, (sy, sx) in enumerate(seeds):
        dy = np.minimum(np.abs(yy - sy), size - np.abs(yy - sy))
        dx = np.minimum(np.abs(xx - sx), size - np.abs(xx - sx))
        d = dx * dx + dy * dy
        mask = d < dist
        lab[mask] = k
        dist[mask] = d[mask]
    return lab


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2027)))
    size = int(p.get("size", 50))
    n_grains = int(p.get("initial_grains", 24))
    mcs = int(p.get("mcs", 30))
    gamma = float(p.get("boundary_energy", 1.0))
    t_eff = float(p.get("effective_temperature", 0.45))
    if size < 8 or n_grains < 2 or n_grains >= size * size or mcs < 1 or gamma <= 0 or t_eff < 0:
        raise ValueError("Invalid Potts grain-growth settings")
    labels = _voronoi_initial(size, n_grains, rng)
    before = labels.copy()
    snapshots = [labels.copy()]
    records = []
    attempts = size * size
    for sweep in range(1, mcs + 1):
        for _ in range(attempts):
            i = int(rng.integers(size)); j = int(rng.integers(size))
            nbrs = periodic_neighbors_2d(i, j, size, size)
            ni, nj = nbrs[int(rng.integers(4))]
            old = int(labels[i, j])
            new = int(labels[ni, nj])
            if old == new:
                continue
            e_old = sum(1 for x, y in nbrs if labels[x, y] != old) * gamma
            e_new = sum(1 for x, y in nbrs if labels[x, y] != new) * gamma
            dE = e_new - e_old
            accept = dE <= 0 or (t_eff > 0 and rng.random() < np.exp(-dE / t_eff))
            if accept:
                labels[i, j] = new
        unique = np.unique(labels)
        records.append({"MCS": sweep, "remaining_orientation_labels": len(unique),
                        "mean_area_per_label_sites": labels.size / len(unique),
                        "grain_boundary_fraction": boundary_fraction(labels)})
        if sweep in {max(1, mcs // 4), max(1, mcs // 2), max(1, 3 * mcs // 4), mcs}:
            snapshots.append(labels.copy())
    df = pd.DataFrame(records)
    series = {
        "Remaining labels": {"x": df["MCS"].to_numpy(), "y": df["remaining_orientation_labels"].to_numpy(),
                             "xlabel": "Monte Carlo steps (MCS)", "ylabel": "Remaining labels"},
        "Grain-boundary fraction": {"x": df["MCS"].to_numpy(), "y": df["grain_boundary_fraction"].to_numpy(),
                                    "xlabel": "Monte Carlo steps (MCS)", "ylabel": "Unlike-neighbour bond fraction"},
    }
    metrics = {"initial_labels": int(len(np.unique(before))),
               "final_labels": int(len(np.unique(labels))),
               "final_mean_area_per_label_sites": float(labels.size / len(np.unique(labels))),
               "final_grain_boundary_fraction": boundary_fraction(labels),
               "total_attempted_site_updates": int(mcs * attempts)}
    return make_result(
        name="grain_growth", title="Grain growth · 2D Potts Monte Carlo",
        method="Metropolis Monte Carlo · q-state Potts model",
        description="Simulates curvature-driven coarsening of a simplified 2D grain-orientation field.",
        before=before, after=labels, data=df, series=series, metrics=metrics,
        visual_kind="labels", before_label="Initial grain map", after_label="Coarsened grain map",
        snapshots=snapshots,
        notes=[
            "The Potts Hamiltonian penalizes unlike nearest-neighbour grain labels; Metropolis updates approximate curvature-driven boundary motion.",
            "MCS is an algorithmic time unit. A physical seconds conversion requires a calibrated grain-boundary mobility and model-specific mapping.",
            "The model is 2D, periodic and isotropic; it does not include anisotropic boundary energy, triple-junction physics, solute drag, pinning or elastic effects.",
        ],
        warnings=["Do not report MCS as seconds or claim quantitative grain-growth kinetics without calibration against measured mobility / grain-size evolution."],
        metadata={"size": size, "initial_grains": n_grains, "mcs": mcs,
                  "boundary_energy_model_units": gamma, "effective_temperature_model_units": t_eff,
                  "seed": int(p.get("seed", 2027))},
    )
