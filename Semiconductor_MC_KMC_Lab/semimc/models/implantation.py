"""Stochastic range/straggle implant profile model.

This is an educational Monte Carlo surrogate, NOT a full binary-collision
approximation (BCA) implementation. It samples calibrated projected ranges
and lateral straggles. The API is intentionally suitable for replacement by
SRIM/TRIM or a proper screened-Coulomb BCA engine.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result

# Illustrative coefficients only; users must calibrate against SRIM/TRIM or SIMS.
_RANGE_COEFF = {"B": 2.10, "P": 1.50, "As": 1.00}
_RANGE_POWER = 1.12


def simulate(p: dict) :
    rng = np.random.default_rng(int(p.get("seed", 1234)))
    species = str(p.get("species", "B"))
    if species not in _RANGE_COEFF:
        raise ValueError("species must be one of B, P, As")
    n_ions = int(p.get("n_ions", 5000))
    energy = float(p.get("energy_keV", 20.0))
    dose = float(p.get("dose_cm2", 1e15))
    tilt = float(p.get("tilt_deg", 7.0))
    thickness = float(p.get("thickness_nm", 300.0))
    bins = int(p.get("bins", 80))
    if n_ions < 100 or energy <= 0 or dose <= 0 or thickness <= 0 or bins < 10:
        raise ValueError("Use n_ions >= 100, positive energy/dose/thickness, and bins >= 10")

    mean_range = _RANGE_COEFF[species] * energy ** _RANGE_POWER
    sigma_range = max(0.18 * mean_range, 1.0)
    z = rng.normal(mean_range, sigma_range, n_ions)
    z = np.clip(z, 0.0, None)
    lateral_sigma = max(0.22 * mean_range, 1.0)
    lateral = z * np.tan(np.deg2rad(tilt)) + rng.normal(0.0, lateral_sigma, n_ions)

    z_edges = np.linspace(0.0, thickness, bins + 1)
    x_extent = max(3.5 * lateral_sigma + thickness * abs(np.tan(np.deg2rad(tilt))), 10.0)
    x_edges = np.linspace(-x_extent, x_extent, 81)
    hist, z_edges_out, x_edges_out = np.histogram2d(z, lateral, bins=(z_edges, x_edges))
    before = np.zeros_like(hist)
    after = hist  # rows are depth bins; columns are lateral bins

    # Implant dose represented in atoms/cm^2 per depth bin; convert nm to cm.
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    counts, _ = np.histogram(z, bins=z_edges)
    dz_cm = np.diff(z_edges) * 1e-7
    concentration = dose * (counts / n_ions) / dz_cm
    data = pd.DataFrame({
        "depth_nm": z_centers,
        "ion_count": counts,
        "concentration_cm-3": concentration,
    })

    # Kinchin-Pease-like damage proxy. It is intentionally not a defect-survival model.
    nuclear_fraction_assumed = float(p.get("nuclear_energy_fraction", 0.12))
    displacement_energy_eV = float(p.get("displacement_energy_eV", 25.0))
    if not (0 < nuclear_fraction_assumed <= 1) or displacement_energy_eV <= 0:
        raise ValueError("Invalid damage proxy parameters")
    damage_mean = nuclear_fraction_assumed * energy * 1000.0 / (2.0 * displacement_energy_eV)
    damage_events = rng.poisson(damage_mean, n_ions)

    profile_mean = float(np.mean(z))
    profile_std = float(np.std(z, ddof=1))
    in_target = float(np.mean(z <= thickness))
    metrics = {
        "sampled_mean_range_nm": profile_mean,
        "sampled_range_straggle_nm": profile_std,
        "fraction_stopped_within_target": in_target,
        "estimated_peak_concentration_cm-3": float(concentration.max(initial=0.0)),
        "mean_damage_proxy_events_per_ion": float(np.mean(damage_events)),
        "nominal_model_range_nm": mean_range,
    }
    data_ions = pd.DataFrame({"ion_id": np.arange(n_ions), "depth_nm": z,
                              "lateral_nm": lateral, "damage_proxy_events": damage_events})
    series = {
        "Implanted concentration profile": {
            "x": z_centers, "y": concentration, "xlabel": "Depth (nm)",
            "ylabel": "Concentration (cm⁻³)"},
        "Ion depth histogram": {
            "x": z_centers, "y": counts, "xlabel": "Depth (nm)", "ylabel": "Ion count"},
    }
    return make_result(
        name="ion_implantation", title="Ion implantation in silicon",
        method="Monte Carlo · stochastic range/straggle surrogate",
        description="Samples ion projected ranges and lateral straggles for an illustrative implant profile.",
        before=before, after=after, data=data, series=series, metrics=metrics,
        visual_kind="implant", before_label="Undoped target (initial)",
        after_label="Sampled implanted-ion density (final)",
        notes=[
            "This module is not a full BCA/TRIM implementation: it samples a Gaussian-like range distribution using explicit illustrative coefficients.",
            "The damage number is a Kinchin–Pease-style proxy from an assumed nuclear energy fraction; it does not predict surviving defects or amorphization.",
            "Dose is normalized into a depth concentration histogram; finite lateral extent is not a full 3D wafer model.",
        ],
        warnings=["Calibrate range coefficients and straggle against SRIM/TRIM and/or experimental SIMS before making quantitative claims."],
        metadata={"species": species, "energy_keV": energy, "dose_cm2": dose,
                  "tilt_deg": tilt, "n_ions": n_ions, "seed": int(p.get("seed", 1234)),
                  "individual_ions": data_ions},
    )
