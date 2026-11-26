# Scientific scope, governing equations and validation plan

This file distinguishes the model actually implemented from a more complete research implementation. A module name or a plausible-looking image does not imply predictive accuracy.

## 1. Ion implantation — Monte Carlo

**Implemented:** A stochastic projected-range / lateral-straggle sampler. Species-dependent illustrative range coefficients generate sampled ion endpoints and histograms. A damage proxy is estimated from an assumed nuclear-energy fraction.

**Not implemented:** Full binary-collision approximation (BCA), screened ZBL trajectory integration, crystal channeling, thermal vibrations, dynamic damage accumulation, electronic/nuclear stopping tables, recoil cascades and amorphization.

For a full BCA model, calculate impact parameters, solve screened two-body scattering, partition energy into nuclear/electronic losses and follow recoils. A common universal screened Coulomb potential uses a ZBL-type screening function. Benchmark against SRIM/TRIM or a crystal-aware BCA code and SIMS data.

Useful observables:

- Projected range: `Rp = mean(z_i)`.
- Range straggling: `sigma_R = std(z_i)`.
- Dose-normalized profile: `C_j = dose * (N_j/N_total) / delta_z_cm`.

The current coefficients are placeholders. Replace them with calibrated species/energy/orientation tables or a real BCA engine.

## 2. Dopant diffusion — Kinetic Monte Carlo

**Implemented:** One mobile vacancy and one dopant tracer on a periodic diamond-cubic network with four nearest neighbours per site. Vacancy migration and vacancy–dopant exchange have separate Arrhenius rates. Residence-time KMC selects events in proportion to rates and advances physical time.

For event i:

`k_i = nu_i * exp(-E_i / (k_B T))`

`K = sum_i(k_i)`

`P_i = k_i/K`

`delta_t = -ln(u)/K`, with `u` uniform on (0,1).

The 3D tracer estimate uses the Einstein relation:

`D = slope(MSD versus t) / 6`.

**Not implemented:** Multiple defect species, charge states, Fermi-level dependence, vacancy/interstitial-mediated mechanisms in parallel, dopant-defect binding/dissociation, clustering, concentration dependence, strain, electric fields or a realistic equilibrium defect population.

The one-vacancy tracer calculation is a useful KMC algorithm demonstration, not a calibrated B/P/As diffusivity prediction. Source migration, exchange, binding and prefactor parameters from dopant- and charge-state-specific literature/first-principles data.

Validation:

1. Verify diamond network has four unique, reciprocal nearest neighbours.
2. Verify KMC time is monotonic and positive.
3. Verify random-walk MSD scales linearly with time in the long-time limit.
4. Compare D(T) to an appropriate analytical/experimental tracer-diffusion reference after adding the correct defect populations and correlation factors.
5. Perform replicate and finite-size convergence studies.

## 3. Thin-film deposition — KMC solid-on-solid model

**Implemented:** Periodic height lattice with deposition, surface hopping and desorption events. Every event has a rate; the total rate determines physical time advancement. The model reports average thickness and RMS roughness.

`k_i = nu_i exp(-E_i/(k_B T))`

`W = sqrt(<h^2> - <h>^2)`

**Not implemented:** Material-specific precursor reaction networks, gas-phase transport, sticking coefficients, surface reconstruction, overhangs, epitaxial registry, plasma species, elastic stress and full PVD/CVD/ALD chemistry.

Default `effective_attempt_frequency_s` is coarse grained. Calibrate flux, sticking, diffusion/desorption barriers, prefactors and bond model for a named process.

Validation: check atom accounting (deposited - desorbed = change in film atoms), positive KMC time, mesh-size convergence, growth-rate comparison and AFM roughness comparison.

## 4. Grain growth — Potts Monte Carlo

**Implemented:** q-state nearest-neighbour Potts model. The Hamiltonian penalizes unlike neighbouring labels:

`H = J * sum_<ij> (1 - delta(s_i,s_j))`.

A trial update is accepted if `delta_H <= 0`, otherwise with probability `exp(-delta_H/T_eff)`.

**Not implemented:** Physical grain-boundary mobility, anisotropic energy, triple-junction mobility, solute drag, Zener pinning, stored-energy driving forces, recrystallization nucleation, 3D topology or calibrated seconds.

One MCS is N attempted updates for N sites. MCS must not be called physical time without a mobility/time calibration. Validate scaling of mean grain size with MCS and compare against known curvature-driven coarsening trends and measured grain-size distributions.

## 5. Binary alloys — Kawasaki Monte Carlo

**Implemented:** Nearest-neighbour atom swaps conserve composition. The Ising pair Hamiltonian is `H = -J sum_<ij> s_i s_j`. Negative J promotes unlike-neighbour ordering on the square lattice; positive J promotes like-neighbour clustering/phase separation.

**Not implemented:** A real alloy's crystal structure, multicomponent interactions, cluster expansion, elastic contributions, vacancies, chemical potentials, phase diagram calibration or physical kinetics. MCS is algorithmic time.

Validate exact A/B composition conservation, energy change against a direct recomputation test, order parameter behaviour and comparison with an exact/known square-lattice model where applicable.

## 6. Defect annealing — object KMC demonstrator

**Implemented:** Vacancy and interstitial objects hop on a periodic 2D lattice and recombine when nearest neighbours. Migration and recombination rates use Arrhenius forms.

**Not implemented:** Crystallographic silicon topology, defect charge state, cluster size distributions, interstitial loops, extended defects, dopant trapping, dissociation networks or realistic reaction barriers. This is not the same as a full semiconductor process OKMC model.

Validate event-rate selection, monotonic time, defect inventory accounting and limiting cases (e.g. no recombination when the recombination rate is set to zero).

## 7. Plasma etching — MC surface-yield prototype

**Implemented:** MC ion starting positions and incidence angles, simple ray shadowing against a mask, and a stochastic effective removal depth.

**Not implemented:** Plasma discharge, sheath and ion energy distributions, neutral/radical transport, ion reflection/scattering, charging, chemical surface reaction network, passivation, redeposition, aspect-ratio-dependent transport, or true MC/KMC coupling.

Validate mask protection, convergence versus ion count, angular trend and calibrated etch rate / SEM profile.

## 8. Thermal oxidation — continuum baseline

**Implemented:** Classical Deal–Grove linear-parabolic equation:

`x_ox^2 + A*x_ox = B*(t + tau)`.

This is not an MC/KMC model. The constants depend on temperature, ambient, pressure, orientation and other process conditions. Validate thickness vs time against oxidation measurements / literature data.

## 9. Photolithography — MC photon absorption demonstrator

**Implemented:** Photon x-position follows an ideal periodic binary mask; absorption depth follows Beer–Lambert exponential sampling. A saturating conversion surrogate generates a thresholded latent image.

**Not implemented:** Vector optical diffraction, partial coherence, standing waves, photoacid generation/diffusion, chemical amplification, post-exposure bake, developer transport, stochastic dissolution and full resist chemistry.

Validate photon statistics against Beer–Lambert attenuation and then extend with calibrated optical and chemical resist models before comparing critical dimensions or line-edge roughness.

## Cross-cutting validation requirements

For every research model, record:

- Material/process identity, crystal orientation and model dimensionality.
- Parameter value, unit, temperature, source/DOI, uncertainty and calibration status.
- Random seed, software version, event count, mesh size and boundary condition.
- Independent replicate statistics and confidence intervals.
- Conservation checks and limiting-case tests.
- Finite-size / time-step / event-count convergence.
- Experimental or analytical reference dataset and a defined error metric.

Suggested metrics include RMSE and normalized RMSE. Avoid MAPE where the reference is zero or near zero. For stochastic results, show uncertainty bands rather than a single curve.
