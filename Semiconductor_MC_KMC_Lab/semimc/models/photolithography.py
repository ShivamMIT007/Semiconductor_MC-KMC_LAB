"""Monte Carlo photon absorption / latent-image demonstrator for photoresist.

The model samples photon absorption depth with Beer-Lambert statistics and
an idealized mask pattern. It does not solve full vector diffraction, resist
chemistry, acid-base kinetics, development or stochastic molecular reactions.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result


def simulate(p: dict):
    rng = np.random.default_rng(int(p.get("seed", 2031)))
    n_photons = int(p.get("photons", 80000))
    nx = int(p.get("pixels_x", 96))
    nz = int(p.get("pixels_z", 48))
    thickness = float(p.get("resist_thickness_um", 1.0))
    alpha = float(p.get("absorption_coefficient_um-1", 2.0))
    pitch = float(p.get("pattern_pitch_um", 0.8))
    line_fraction = float(p.get("open_fraction", 0.5))
    threshold = float(p.get("conversion_threshold", 0.35))
    photon_yield = float(p.get("conversion_per_absorbed_photon", 0.02))
    if n_photons < 100 or nx < 16 or nz < 8 or min(thickness, alpha, pitch) <= 0 or not 0 < line_fraction < 1 or not 0 < threshold < 1 or photon_yield <= 0:
        raise ValueError("Invalid photolithography parameters")
    x = rng.uniform(0.0, nx * pitch / 16.0, n_photons)
    # Ideal periodic binary mask: a fixed fraction of each pitch is transparent.
    local = np.mod(x, pitch)
    is_open = local < pitch * line_fraction
    x = x[is_open]
    # Beer-Lambert absorption depth: P(z) = alpha exp(-alpha z).
    u = np.maximum(rng.random(x.size), np.finfo(float).tiny)
    z = -np.log(u) / alpha
    in_resist = z <= thickness
    x = x[in_resist]
    z = z[in_resist]
    x_edges = np.linspace(0.0, nx * pitch / 16.0, nx + 1)
    z_edges = np.linspace(0.0, thickness, nz + 1)
    counts, _, _ = np.histogram2d(z, x, bins=(z_edges, x_edges))
    # Poisson-like molecular conversion surrogate per voxel.
    conversion = 1.0 - np.exp(-photon_yield * counts)
    latent = (conversion >= threshold).astype(float)
    before = np.zeros_like(latent)
    after = latent.copy()
    x_centres = 0.5 * (x_edges[:-1] + x_edges[1:])
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    rows = []
    for iz, zc in enumerate(z_centres):
        for ix, xc in enumerate(x_centres):
            rows.append((xc, zc, counts[iz, ix], conversion[iz, ix], latent[iz, ix]))
    df = pd.DataFrame(rows, columns=["x_um", "depth_um", "absorbed_photons", "conversion_fraction_surrogate", "above_threshold"])
    metrics = {"incident_photons": n_photons,
               "absorbed_photons_in_resist": int(x.size),
               "fraction_absorbed_in_resist": float(x.size / n_photons),
               "mean_conversion_surrogate": float(conversion.mean()),
               "fraction_voxels_above_threshold": float(latent.mean()),
               "resist_thickness_um": thickness}
    profile = conversion.mean(axis=1)
    series = {"Mean conversion versus depth": {"x": z_centres, "y": profile,
                                                "xlabel": "Resist depth (µm)", "ylabel": "Mean conversion surrogate"}}
    return make_result(
        name="photolithography", title="Photolithography · Monte Carlo photon absorption",
        method="Monte Carlo photon absorption · Beer–Lambert sampling",
        description="Samples photon absorption depths beneath an ideal periodic binary mask to form a stochastic latent image.",
        before=before, after=after, data=df, series=series, metrics=metrics,
        visual_kind="implant", before_label="Unexposed resist", after_label="Thresholded latent-image surrogate",
        notes=[
            "Absorption depth follows an exponential Beer–Lambert distribution and photon x-positions follow an ideal binary mask.",
            "The conversion fraction is an illustrative saturating function of absorbed photon count, not a calibrated chemically amplified resist model.",
            "Optical diffraction, standing waves, acid diffusion, post-exposure bake chemistry and developer transport are omitted.",
        ],
        warnings=["Use a calibrated optical/resist model and SEM/CD data before treating predicted critical dimensions as quantitative."],
        metadata={"photons": n_photons, "absorption_coefficient_um-1": alpha,
                  "pattern_pitch_um": pitch, "open_fraction": line_fraction,
                  "conversion_threshold": threshold, "seed": int(p.get("seed", 2031))},
    )
