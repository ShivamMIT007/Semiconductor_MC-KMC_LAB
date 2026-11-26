# Semiconductor Monte Carlo / Kinetic Monte Carlo Lab

A VS Code + Python + Streamlit workbench containing separate simulation modules for semiconductor processing and materials microstructure. Each run can produce initial/final state images, time-series graphs, a CSV dataset, an animation (where snapshots are available), numerical metrics, metadata and a ZIP research record.

> **Scope warning:** This is a research-oriented *prototype and teaching framework*, not a validated commercial process simulator. Some modules are deliberately simplified. Read `SCIENTIFIC_SCOPE.md` before interpreting results quantitatively. All parameters labelled placeholder/effective/illustrative must be sourced and calibrated for a named material and process before publication or process decisions.

## Included models

| Process | Method implemented | Main outputs | Important scope note |
|---|---|---|---|
| Ion implantation in Si | MC range/straggle surrogate | Sampled depth/lateral profiles, range, straggle, damage proxy | Not a full screened-Coulomb BCA / SRIM replacement |
| Vacancy-mediated dopant tracer diffusion | Residence-time KMC on explicit 3D diamond-cubic network | Physical time, tracer MSD, estimated tracer D | One vacancy + one tracer per run; default barriers are placeholders |
| Thin-film deposition | KMC solid-on-solid surface growth | Height maps, thickness, roughness, event counts | Coarse-grained model; not chemically specific PVD/CVD/ALD |
| Grain growth | Metropolis q-state Potts MC | Before/after grain maps, boundary fraction, orientation-label count | MCS is algorithmic time, not seconds |
| Binary alloy ordering | Conserved Kawasaki MC, antiferromagnetic Ising coupling | Ordering parameter, energy, composition | Generic square lattice, not a fitted alloy Hamiltonian |
| Binary alloy phase separation | Conserved Kawasaki MC, ferromagnetic Ising coupling | Cluster fraction, energy, composition | Generic square lattice, not a real alloy phase diagram |
| Vacancy/interstitial annealing | Object KMC on a coarse 2D lattice | Defect populations, recombination count, physical event time | Demonstrator; no full Si defect species / charge states / clusters |
| Plasma etching | MC ion incidence + effective stochastic removal | Etch-depth profile, hit fraction, nonuniformity | Not a plasma chemistry, sheath or full hybrid MC-KMC solver |
| Thermal oxidation | Deal-Grove continuum equation | Oxide thickness vs time | Correctly classified as continuum, not MC/KMC |
| Photolithography | MC Beer-Lambert photon absorption | Exposure / conversion map, depth profile | Does not solve diffraction or full resist chemistry / development |

## Windows / VS Code installation

1. Download and extract the project ZIP.
2. In VS Code choose **File → Open Folder** and open the folder that directly contains `app.py`, `requirements.txt`, and `semimc/`.
3. Open **Terminal → New Terminal**.
4. Confirm the current folder contains the files:

   ```powershell
   dir
   Test-Path .\requirements.txt
   ```

   `Test-Path` must return `True`. If it returns `False`, you are in the outer extraction folder; `cd` into the inner `Semiconductor_MC_KMC_Lab` folder.

5. Create an environment and install packages:

   ```powershell
   py -m venv .venv
   .\.venv\Scripts\python.exe -m pip install --upgrade pip
   .\.venv\Scripts\python.exe -m pip install -r requirements.txt
   ```

6. Launch the website:

   ```powershell
   .\.venv\Scripts\python.exe -m streamlit run app.py
   ```

7. Open the local URL printed by Streamlit, normally `http://localhost:8501`.

Alternatively, double-click `run_app.bat` from File Explorer. First installation needs internet access to download Python packages.

## Run models directly from Python

Example:

```python
from semimc.models.dopant_diffusion import simulate

result = simulate({
    "n_cells": 3,
    "temperature_K": 1100.0,
    "vacancy_migration_eV": 0.50,
    "dopant_exchange_eV": 0.55,
    "attempt_frequency_s": 1.0e13,
    "events_per_run": 5000,
    "replicates": 8,
    "seed": 2025,
})
print(result.metrics)
print(result.data.head())
```

Model entry points and parameters are registered in `semimc/catalog.py`. Individual algorithms live under `semimc/models/`.

## Generate example outputs

```powershell
.\.venv\Scripts\python.exe scripts\generate_examples.py
```

The script writes a folder per model under `examples/`, including `summary.json`, `dataset.csv`, before/after PNGs, a graph PNG and an evolution GIF when the model stores snapshots.

## Run automated checks

Install pytest if desired:

```powershell
.\.venv\Scripts\python.exe -m pip install pytest
.\.venv\Scripts\python.exe -m pytest -q
```

Or run the test file with Python's standard library:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

Tests verify reproducibility, lattice coordination/reciprocity, positive finite outputs, alloy composition conservation, defect-count behaviour and monotonic oxidation thickness. They are code checks, not experimental validation.

## Output and export

The app provides:

- Initial and final images.
- Time-series graphs.
- CSV numerical datasets.
- Event-driven animation GIFs if snapshots exist.
- Numerical metrics and model metadata.
- ZIP run bundle with CSV, PNG, GIF (when available), and JSON.

Runs can also be saved to the local `results/` directory from the Export tab.

## Suggested research sequence

1. Benchmark the implantation sampler against SRIM/TRIM and SIMS; replace the current range surrogate with BCA before making process predictions.
2. Extend the diamond-lattice vacancy tracer into a calibrated object-KMC network with V, I, dopant-defect pairs, charge states, binding and dissociation events.
3. Couple the as-implanted defect distribution to the annealing/diffusion KMC initial condition.
4. Calibrate surface KMC event networks for a named PVD/CVD/ALD chemistry.
5. Extend plasma etching with ion energy/angular distributions, neutral/radical transport, surface KMC chemistry, charging and profile evolution.
6. Compare each model with experiment and report uncertainty, convergence and parameter provenance.

See `SCIENTIFIC_SCOPE.md` and `REFERENCES.md` for model limitations and literature starting points.
