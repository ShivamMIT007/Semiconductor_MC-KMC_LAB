# Dataset and unit conventions

Each web-app run exports a CSV table. The exact columns depend on the selected model.

| Model | Principal columns | Units / interpretation |
|---|---|---|
| Ion implantation | `depth_nm`, `ion_count`, `concentration_cm-3` | Depth in nm; concentration derived from user dose and sampled histogram |
| Dopant diffusion | `normalized_time`, `mean_physical_time_s`, `ensemble_mean_tracer_MSD_nm2` | Physical KMC time in s; MSD in nm² |
| Thin-film deposition | `event`, `time_s`, `mean_thickness_ML`, `surface_roughness_ML` | KMC event index, physical time, monolayers and RMS monolayer height |
| Grain growth | `MCS`, `remaining_orientation_labels`, `grain_boundary_fraction` | Algorithmic Potts sweeps; boundary fraction is dimensionless |
| Binary alloy | `MCS`, `energy_eV_per_site`, `unlike_neighbour_fraction` | Algorithmic Kawasaki sweeps; energy uses user interaction parameter |
| Defect annealing | `event`, `time_s`, `vacancy_count`, `interstitial_count` | Physical KMC time in s; object counts |
| Plasma etching | `lateral_position_nm`, `etch_depth_nm` | Model lateral coordinate and effective depth; lateral pitch mapping is not built in |
| Thermal oxidation | `time_s`, `oxide_thickness_nm` | Seconds and nm |
| Photolithography | `x_um`, `depth_um`, `absorbed_photons`, `conversion_fraction_surrogate` | Micrometres, photon counts and uncalibrated conversion surrogate |

Do not compare values across models unless their units, normalization, boundary conditions and parameter definitions match.
