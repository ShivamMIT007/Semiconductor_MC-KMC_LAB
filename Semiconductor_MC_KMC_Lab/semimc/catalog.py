"""UI catalogue: model implementations, default parameters and control metadata."""
from semimc.models.implantation import simulate as implantation
from semimc.models.dopant_diffusion import simulate as dopant_diffusion
from semimc.models.thin_film import simulate as thin_film
from semimc.models.grain_growth import simulate as grain_growth
from semimc.models.alloy import simulate_ordering, simulate_phase_separation
from semimc.models.defect_annealing import simulate as defect_annealing
from semimc.models.plasma_etch import simulate as plasma_etch
from semimc.models.oxidation import simulate as oxidation
from semimc.models.photolithography import simulate as photolithography

# (key, label, widget type, default, min, max, step, help)
CATALOG = {
    "Ion implantation (MC)": {
        "fn": implantation, "params": [
            ("species", "Implant species", "select", "B", None, None, None, "Boron, phosphorus or arsenic."),
            ("energy_keV", "Ion energy (keV)", "float", 20.0, 1.0, 500.0, 1.0, "Incident ion energy."),
            ("dose_cm2", "Dose (ions/cm²)", "float", 1e15, 1e10, 1e17, 1e13, "Implanted fluence."),
            ("n_ions", "Histories / ions", "int", 5000, 100, 100000, 100, "Number of sampled ion histories."),
            ("tilt_deg", "Tilt angle (degrees)", "float", 7.0, 0.0, 45.0, 0.5, "Mean implant tilt away from surface normal."),
            ("thickness_nm", "Target thickness (nm)", "float", 300.0, 20.0, 2000.0, 10.0, "Depth range represented in the histogram."),
            ("bins", "Depth bins", "int", 80, 10, 300, 5, "Histogram resolution."),
            ("seed", "Random seed", "int", 1234, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Dopant diffusion (KMC)": {
        "fn": dopant_diffusion, "params": [
            ("n_cells", "Diamond unit cells per direction", "int", 3, 2, 8, 1, "Periodic diamond-cubic network size; sites = 8 n³."),
            ("temperature_K", "Temperature (K)", "float", 1100.0, 300.0, 1800.0, 10.0, "Absolute temperature."),
            ("vacancy_migration_eV", "Vacancy migration barrier (eV)", "float", 0.50, 0.01, 3.0, 0.01, "Placeholder barrier; replace with a source for the specific defect."),
            ("dopant_exchange_eV", "Dopant-vacancy exchange barrier (eV)", "float", 0.55, 0.01, 3.0, 0.01, "Placeholder exchange barrier."),
            ("attempt_frequency_s", "Attempt frequency (s⁻¹)", "float", 1e13, 1e6, 1e14, 1e11, "Harmonic prefactor / effective event frequency."),
            ("events_per_run", "KMC events per replicate", "int", 5000, 100, 100000, 100, "Number of vacancy events in each run."),
            ("replicates", "Independent replicates", "int", 8, 2, 50, 1, "Used for a stochastic uncertainty estimate."),
            ("seed", "Random seed", "int", 2025, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Thin-film deposition (KMC)": {
        "fn": thin_film, "params": [
            ("size", "Surface size (sites per side)", "int", 18, 6, 40, 1, "Periodic square surface."),
            ("events", "KMC events", "int", 6000, 100, 30000, 100, "Total attempted KMC events."),
            ("temperature_K", "Temperature (K)", "float", 650.0, 200.0, 1500.0, 10.0, "Absolute temperature."),
            ("flux_ml_s", "Deposition flux (ML/s/site)", "float", 0.25, 0.001, 20.0, 0.01, "Arrival rate at each surface column."),
            ("effective_attempt_frequency_s", "Effective attempt frequency (s⁻¹)", "float", 1e5, 1e1, 1e13, 1e3, "Coarse-grained default; calibrate for atomistic use."),
            ("diffusion_barrier_eV", "Surface diffusion barrier (eV)", "float", 0.65, 0.01, 3.0, 0.01, "Activation barrier for hopping."),
            ("bond_penalty_eV", "Lateral bond penalty (eV)", "float", 0.08, 0.0, 2.0, 0.01, "Extra barrier per lateral bond."),
            ("desorption_barrier_eV", "Desorption barrier (eV)", "float", 1.8, 0.01, 5.0, 0.01, "Activation barrier for desorption."),
            ("seed", "Random seed", "int", 2026, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Grain growth (Potts MC)": {
        "fn": grain_growth, "params": [
            ("size", "Grid size", "int", 50, 10, 120, 2, "2D periodic grain-orientation grid."),
            ("initial_grains", "Initial grain seeds", "int", 24, 2, 100, 1, "Initial Voronoi grain count."),
            ("mcs", "Monte Carlo steps (MCS)", "int", 30, 1, 300, 1, "One MCS = one attempted update per lattice site."),
            ("boundary_energy", "Boundary energy (model units)", "float", 1.0, 0.01, 10.0, 0.05, "Potts unlike-neighbour penalty."),
            ("effective_temperature", "Effective temperature (model units)", "float", 0.45, 0.0, 5.0, 0.05, "Metropolis noise parameter; not Kelvin."),
            ("seed", "Random seed", "int", 2027, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Binary alloy ordering (Kawasaki MC)": {
        "fn": simulate_ordering, "params": [
            ("size", "Grid size", "int", 40, 8, 100, 2, "2D square lattice."),
            ("mcs", "Monte Carlo steps (MCS)", "int", 30, 1, 200, 1, "Conserved nearest-neighbour exchanges."),
            ("temperature_K", "Temperature (K)", "float", 180.0, 50.0, 2000.0, 10.0, "Absolute temperature."),
            ("interaction_eV", "Interaction J (eV)", "float", -0.045, -0.5, 0.5, 0.005, "Negative J favours unlike neighbours in E=-JΣsisj."),
            ("fraction_A", "Initial A fraction", "float", 0.5, 0.05, 0.95, 0.01, "Composition is conserved by Kawasaki swaps."),
            ("seed", "Random seed", "int", 2028, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Binary alloy phase separation (Kawasaki MC)": {
        "fn": simulate_phase_separation, "params": [
            ("size", "Grid size", "int", 40, 8, 100, 2, "2D square lattice."),
            ("mcs", "Monte Carlo steps (MCS)", "int", 30, 1, 200, 1, "Conserved nearest-neighbour exchanges."),
            ("temperature_K", "Temperature (K)", "float", 300.0, 50.0, 2000.0, 10.0, "Absolute temperature."),
            ("interaction_eV", "Interaction J (eV)", "float", 0.05, 0.001, 0.5, 0.005, "Positive J favours like-neighbour clustering."),
            ("fraction_A", "Initial A fraction", "float", 0.5, 0.05, 0.95, 0.01, "Composition is conserved by Kawasaki swaps."),
            ("seed", "Random seed", "int", 2028, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Defect annealing (object KMC)": {
        "fn": defect_annealing, "params": [
            ("size", "Grid size", "int", 36, 8, 100, 2, "2D periodic coarse defect lattice."),
            ("initial_pairs", "Initial vacancy-interstitial pairs", "int", 45, 1, 500, 1, "Equal starting populations."),
            ("events", "KMC event limit", "int", 5000, 100, 30000, 100, "Maximum event count."),
            ("temperature_K", "Temperature (K)", "float", 900.0, 200.0, 1800.0, 10.0, "Absolute temperature."),
            ("vacancy_migration_eV", "Vacancy migration barrier (eV)", "float", 0.65, 0.01, 3.0, 0.01, "Placeholder parameter."),
            ("interstitial_migration_eV", "Interstitial migration barrier (eV)", "float", 0.25, 0.01, 3.0, 0.01, "Placeholder parameter."),
            ("recombination_barrier_eV", "Recombination barrier (eV)", "float", 0.05, 0.0, 3.0, 0.01, "Placeholder parameter."),
            ("attempt_frequency_s", "Migration prefactor (s⁻¹)", "float", 1e12, 1e6, 1e14, 1e10, "Event-rate prefactor."),
            ("recombination_prefactor_s", "Recombination prefactor (s⁻¹)", "float", 1e12, 1e6, 1e14, 1e10, "Recombination-rate prefactor."),
            ("seed", "Random seed", "int", 2029, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Plasma etching (MC surface-yield prototype)": {
        "fn": plasma_etch, "params": [
            ("columns", "Lateral columns", "int", 100, 20, 300, 5, "2D cross-section resolution."),
            ("ions", "Incident ion histories", "int", 20000, 100, 100000, 100, "Number of sampled ions."),
            ("ion_energy_eV", "Ion energy (eV)", "float", 200.0, 10.0, 2000.0, 10.0, "Representative incident energy."),
            ("angular_sigma_deg", "Angular spread σ (degrees)", "float", 8.0, 0.0, 45.0, 0.5, "Gaussian incidence angle spread."),
            ("removal_nm_per_ion", "Effective removal (nm/ion hit)", "float", 0.002, 0.0, 0.1, 0.0005, "Calibrate to measured etch rate."),
            ("mask_fraction_each_side", "Mask fraction on each side", "float", 0.15, 0.0, 0.4, 0.01, "Fraction of columns protected at each edge."),
            ("seed", "Random seed", "int", 2030, 0, 9999999, 1, "Reproducible random stream."),
        ]},
    "Thermal oxidation (Deal–Grove continuum)": {
        "fn": oxidation, "params": [
            ("A_nm", "Linear-parabolic constant A (nm)", "float", 120.0, 0.1, 5000.0, 1.0, "Ambient/temperature-specific constant."),
            ("B_nm2_s", "Parabolic rate constant B (nm²/s)", "float", 18.0, 0.001, 10000.0, 0.1, "Ambient/temperature-specific constant."),
            ("initial_time_shift_s", "Initial time shift τ (s)", "float", 0.0, 0.0, 1e6, 1.0, "Initial-condition correction."),
            ("duration_s", "Oxidation duration (s)", "float", 3600.0, 1.0, 1e7, 60.0, "Total oxidation time."),
            ("initial_oxide_nm", "Initial oxide (nm)", "float", 0.0, 0.0, 1000.0, 1.0, "Starting oxide thickness."),
            ("points", "Time samples", "int", 150, 10, 1000, 10, "Output resolution."),
        ]},
    "Photolithography (MC photon absorption)": {
        "fn": photolithography, "params": [
            ("photons", "Photon histories", "int", 80000, 100, 500000, 1000, "Monte Carlo photon count."),
            ("pixels_x", "Lateral pixels", "int", 96, 16, 300, 8, "Image resolution in x."),
            ("pixels_z", "Depth pixels", "int", 48, 8, 200, 4, "Image resolution in depth."),
            ("resist_thickness_um", "Resist thickness (µm)", "float", 1.0, 0.1, 10.0, 0.1, "Photoresist thickness."),
            ("absorption_coefficient_um-1", "Absorption coefficient (µm⁻¹)", "float", 2.0, 0.01, 50.0, 0.1, "Beer–Lambert parameter."),
            ("pattern_pitch_um", "Mask pitch (µm)", "float", 0.8, 0.1, 10.0, 0.1, "Ideal repeating mask period."),
            ("open_fraction", "Transparent fraction", "float", 0.5, 0.05, 0.95, 0.01, "Open portion of each mask pitch."),
            ("conversion_threshold", "Latent-image threshold", "float", 0.35, 0.01, 0.99, 0.01, "Threshold for visualization only."),
            ("seed", "Random seed", "int", 2031, 0, 9999999, 1, "Reproducible random stream."),
        ]},
}
