"""Conserved-composition binary alloy lattice Monte Carlo models.

Kawasaki nearest-neighbour swaps conserve A/B composition. Antiferromagnetic
coupling gives a checkerboard ordering tendency; ferromagnetic coupling gives
like-species clustering / phase separation. Energies are model parameters,
not fitted values for a specific named alloy.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from collections import deque
from semimc.common import make_result, boundary_fraction, periodic_neighbors_2d
K_B_EV = 8.617333262145e-5


def _energy(spins: np.ndarray, J: float) -> float:
    # E = -J sum over each right/down bond once.
    return float(-J * (np.sum(spins * np.roll(spins, -1, axis=0)) +
                       np.sum(spins * np.roll(spins, -1, axis=1))))


def _largest_cluster_fraction(spins: np.ndarray) -> float:
    n, m = spins.shape
    seen = np.zeros((n, m), dtype=bool)
    largest = 0
    for i in range(n):
        for j in range(m):
            if seen[i, j]:
                continue
            val = spins[i, j]
            q = deque([(i, j)])
            seen[i, j] = True
            count = 0
            while q:
                x, y = q.popleft(); count += 1
                for xx, yy in periodic_neighbors_2d(x, y, n, m):
                    if not seen[xx, yy] and spins[xx, yy] == val:
                        seen[xx, yy] = True
                        q.append((xx, yy))
            largest = max(largest, count)
    return largest / spins.size


def simulate(p: dict, mode: str = "ordering"):
    rng = np.random.default_rng(int(p.get("seed", 2028)))
    size = int(p.get("size", 40))
    mcs = int(p.get("mcs", 30))
    temp = float(p.get("temperature_K", 180.0 if mode == "ordering" else 300.0))
    J = float(p.get("interaction_eV", -0.045 if mode == "ordering" else 0.05))
    fraction_a = float(p.get("fraction_A", 0.5))
    if size < 6 or mcs < 1 or temp <= 0 or J == 0 or not 0 < fraction_a < 1:
        raise ValueError("Invalid alloy simulation parameters")
    nsites = size * size
    n_a = int(round(fraction_a * nsites))
    spins = np.ones(nsites, dtype=np.int8)
    spins[n_a:] = -1
    rng.shuffle(spins)
    spins = spins.reshape(size, size)
    before = spins.copy()
    checker = np.fromfunction(lambda i, j: np.where((i + j) % 2 == 0, 1, -1), (size, size), dtype=int)
    records = []
    snapshots = [spins.copy()]
    n_attempts = nsites
    beta = 1.0 / (K_B_EV * temp)
    for sweep in range(1, mcs + 1):
        for _ in range(n_attempts):
            i = int(rng.integers(size)); j = int(rng.integers(size))
            ni, nj = periodic_neighbors_2d(i, j, size, size)[int(rng.integers(4))]
            si = int(spins[i, j]); sj = int(spins[ni, nj])
            if si == sj:
                continue
            # Local energy change for swapping unlike neighbours; their mutual bond is unchanged.
            sum_i = sum(int(spins[x, y]) for x, y in periodic_neighbors_2d(i, j, size, size) if (x, y) != (ni, nj))
            sum_j = sum(int(spins[x, y]) for x, y in periodic_neighbors_2d(ni, nj, size, size) if (x, y) != (i, j))
            dE = -J * ((sj - si) * sum_i + (si - sj) * sum_j)
            if dE <= 0 or rng.random() < np.exp(-dE * beta):
                spins[i, j], spins[ni, nj] = sj, si
        energy_per_site = _energy(spins, J) / nsites
        unlike = boundary_fraction(spins)
        if mode == "ordering":
            order = float(abs(np.mean(spins * checker)))
            metric_name = "long_range_checkerboard_order_parameter"
            metric_value = order
        else:
            order = _largest_cluster_fraction(spins)
            metric_name = "largest_same_species_cluster_fraction"
            metric_value = order
        records.append({"MCS": sweep, "energy_eV_per_site": energy_per_site,
                        "unlike_neighbour_fraction": unlike, metric_name: metric_value,
                        "composition_A": float(np.mean(spins == 1))})
        if sweep in {max(1, mcs // 4), max(1, mcs // 2), max(1, 3 * mcs // 4), mcs}:
            snapshots.append(spins.copy())
    df = pd.DataFrame(records)
    order_col = metric_name
    series = {
        "Energy per site": {"x": df["MCS"].to_numpy(), "y": df["energy_eV_per_site"].to_numpy(),
                            "xlabel": "Monte Carlo steps (MCS)", "ylabel": "Energy (eV/site)"},
        "Order / cluster metric": {"x": df["MCS"].to_numpy(), "y": df[order_col].to_numpy(),
                                   "xlabel": "Monte Carlo steps (MCS)", "ylabel": order_col.replace("_", " ")},
    }
    return make_result(
        name=f"binary_alloy_{mode}", title="Binary alloy " + ("ordering" if mode == "ordering" else "phase separation"),
        method="Metropolis Kawasaki Monte Carlo · conserved composition",
        description="Nearest-neighbour A/B atom exchanges under a pair-interaction Ising Hamiltonian.",
        before=before, after=spins.copy(), data=df, series=series,
        metrics={"initial_A_fraction": float(np.mean(before == 1)),
                 "final_A_fraction": float(np.mean(spins == 1)),
                 "final_energy_eV_per_site": float(_energy(spins, J) / nsites),
                 metric_name: metric_value,
                 "final_unlike_neighbour_fraction": boundary_fraction(spins)},
        visual_kind="spins", before_label="Random alloy", after_label="Evolved alloy configuration",
        snapshots=snapshots,
        notes=[
            "Kawasaki exchanges conserve the number of A and B atoms exactly.",
            "For ordering, a negative J in E=-JΣsᵢsⱼ favours unlike neighbours on the square lattice; the checkerboard order parameter is illustrative.",
            "For phase separation, positive J favours like-neighbour clustering. The square-lattice model is not a crystallographically specific alloy Hamiltonian.",
            "MCS is algorithmic time, not seconds.",
        ],
        warnings=["Interaction energy and lattice topology must be calibrated to a real alloy before predicting a material phase diagram or transition temperature."],
        metadata={"mode": mode, "size": size, "mcs": mcs, "temperature_K": temp,
                  "interaction_eV": J, "fraction_A": fraction_a, "seed": int(p.get("seed", 2028))},
    )


def simulate_ordering(p: dict):
    return simulate(p, "ordering")


def simulate_phase_separation(p: dict):
    return simulate(p, "phase_separation")
