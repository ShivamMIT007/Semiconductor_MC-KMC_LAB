"""Monte Carlo ion-incidence / stochastic-yield plasma-etch profile demonstrator.

This is a 2D cross-section toy model: ions are sampled with an angular
spread, rays are shadow-traced against a masked surface, and a calibrated
effective removal depth is applied at the impact column. Gas-phase plasma
chemistry, charging, feature-scale neutral transport and a surface KMC event
network are not implemented.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2030)))
    nx = int(p.get("columns", 100))
    n_ions = int(p.get("ions", 20000))
    energy = float(p.get("ion_energy_eV", 200.0))
    angle_sigma = float(p.get("angular_sigma_deg", 8.0))
    base_removal = float(p.get("removal_nm_per_ion", 0.002))
    e_ref = float(p.get("reference_energy_eV", 200.0))
    exponent = float(p.get("energy_exponent", 0.5))
    angle_power = float(p.get("angular_power", 1.5))
    mask_fraction = float(p.get("mask_fraction_each_side", 0.15))
    if nx < 20 or n_ions < 100 or energy <= 0 or angle_sigma < 0 or base_removal < 0 or e_ref <= 0 or not 0 <= mask_fraction < 0.45:
        raise ValueError("Invalid plasma etch parameters")
    depth = np.zeros(nx, dtype=float)
    initial = depth.copy()
    mask = np.zeros(nx, dtype=bool)
    nmask = int(nx * mask_fraction)
    mask[:nmask] = True
    if nmask:
        mask[-nmask:] = True
    z_start = 20.0
    dz = 0.5
    etched_hits = 0
    shadowed_hits = 0
    # Each ion is traced from above the wafer towards the surface.
    for _ in range(n_ions):
        x0 = rng.uniform(0, nx - 1)
        theta = np.deg2rad(rng.normal(0.0, angle_sigma))
        hit = None
        for z in np.arange(z_start, -float(depth.max()) - 2.0, -dz):
            x = int(round(x0 + (z_start - z) * np.tan(theta)))
            if x < 0 or x >= nx:
                break
            if mask[x] and z <= 0.0:
                shadowed_hits += 1
                hit = "mask"
                break
            if (not mask[x]) and z <= -depth[x]:
                hit = x
                break
        if isinstance(hit, int):
            c = max(float(np.cos(theta)), 0.0)
            mean_remove = base_removal * (energy / e_ref) ** exponent * c ** angle_power
            depth[hit] += rng.exponential(mean_remove) if mean_remove > 0 else 0.0
            etched_hits += 1
    x_nm = np.arange(nx, dtype=float)
    df = pd.DataFrame({"lateral_column": np.arange(nx), "lateral_position_nm": x_nm,
                       "mask": mask, "etch_depth_nm": depth})
    open_depth = depth[~mask]
    avg_depth = float(open_depth.mean()) if open_depth.size else 0.0
    uniformity = float(np.std(open_depth) / avg_depth) if avg_depth > 0 else 0.0
    series = {
        "Etch profile": {"x": x_nm, "y": depth, "xlabel": "Lateral position (model units)", "ylabel": "Etch depth (nm)"},
    }
    metrics = {"mean_open_area_etch_depth_nm": avg_depth,
               "maximum_etch_depth_nm": float(depth.max()),
               "open_area_depth_nonuniformity_coefficient": uniformity,
               "ions_launched": n_ions, "ions_hitting_open_surface": etched_hits,
               "ions_intercepted_by_mask": shadowed_hits,
               "open_area_hit_fraction": float(etched_hits / n_ions)}
    return make_result(
        name="plasma_etch", title="Plasma etching · stochastic ion-incidence profile model",
        method="Monte Carlo particle incidence + stochastic effective surface removal",
        description="Samples ion incidence angles and shadowing against a masked 2D feature profile.",
        before=initial.reshape(1, -1), after=depth.reshape(1, -1), data=df, series=series,
        metrics=metrics, visual_kind="profile", before_label="Initial flat surface", after_label="Etched depth profile",
        notes=[
            "Ion positions and angles are Monte Carlo sampled; the local removal depth is sampled from an effective exponential yield law.",
            "This is not a plasma chemistry solver and not a complete MC-KMC coupling. It omits gas-phase collisions, sheath energy distributions, charging, neutral radicals and passivation.",
            "Lateral coordinates are model columns; set a physical pitch mapping before comparing to SEM dimensions.",
        ],
        warnings=["Removal-per-ion is an effective user parameter, not a universal sputter yield. Calibrate to etch-rate and SEM cross-section data."],
        metadata={"ion_energy_eV": energy, "angular_sigma_deg": angle_sigma,
                  "removal_nm_per_ion": base_removal, "columns": nx, "ions": n_ions,
                  "seed": int(p.get("seed", 2030))},
    )
