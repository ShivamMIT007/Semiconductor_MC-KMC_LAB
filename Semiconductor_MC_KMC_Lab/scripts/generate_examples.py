from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from semimc.catalog import CATALOG
from semimc.visuals import save_result

OUT = ROOT / "examples"
OUT.mkdir(exist_ok=True)

SMALL_RUNS = {
    "Ion implantation (MC)": {"n_ions": 3000, "bins": 75},
    "Dopant diffusion (KMC)": {"events_per_run": 1800, "replicates": 5},
    "Thin-film deposition (KMC)": {"size": 14, "events": 1800},
    "Grain growth (Potts MC)": {"size": 40, "mcs": 20},
    "Binary alloy ordering (Kawasaki MC)": {"size": 32, "mcs": 20},
    "Binary alloy phase separation (Kawasaki MC)": {"size": 32, "mcs": 20},
    "Defect annealing (object KMC)": {"size": 30, "initial_pairs": 35, "events": 1800},
    "Plasma etching (MC surface-yield prototype)": {"columns": 90, "ions": 9000},
    "Thermal oxidation (Deal–Grove continuum)": {"points": 120},
    "Photolithography (MC photon absorption)": {"photons": 45000, "pixels_x": 80, "pixels_z": 40},
}

for name, overrides in SMALL_RUNS.items():
    print(f"Generating {name} ...")
    params = {item[0]: item[3] for item in CATALOG[name]["params"]}
    params.update(overrides)
    result = CATALOG[name]["fn"](params)
    folder = save_result(result, str(OUT))
    print(f"  saved {folder}")
