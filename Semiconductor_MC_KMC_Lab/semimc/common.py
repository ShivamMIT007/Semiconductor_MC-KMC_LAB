from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any
import numpy as np
import pandas as pd


@dataclass
class SimulationResult:
    name: str
    title: str
    method: str
    description: str
    before: np.ndarray
    after: np.ndarray
    data: pd.DataFrame
    series: dict[str, dict[str, Any]]
    metrics: dict[str, Any]
    notes: list[str]
    warnings: list[str]
    before_label: str = "Initial state"
    after_label: str = "Final state"
    visual_kind: str = "lattice"
    snapshots: list[np.ndarray] | None = None
    metadata: dict[str, Any] | None = None

    def summary(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "title": self.title,
            "method": self.method,
            "description": self.description,
            "metrics": _json_safe(self.metrics),
            "notes": self.notes,
            "warnings": self.warnings,
            "metadata": _json_safe(self.metadata or {}),
        }


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    return str(obj)


def make_result(**kwargs: Any) -> SimulationResult:
    """Construct a result with normalized arrays and default metadata."""
    kwargs.setdefault("warnings", [])
    kwargs.setdefault("notes", [])
    kwargs.setdefault("series", {})
    kwargs.setdefault("metrics", {})
    kwargs.setdefault("metadata", {})
    kwargs.setdefault("snapshots", None)
    if "data" not in kwargs:
        kwargs["data"] = pd.DataFrame()
    kwargs["before"] = np.asarray(kwargs["before"])
    kwargs["after"] = np.asarray(kwargs["after"])
    return SimulationResult(**kwargs)


def periodic_neighbors_2d(i: int, j: int, nx: int, ny: int) -> list[tuple[int, int]]:
    return [((i - 1) % nx, j), ((i + 1) % nx, j), (i, (j - 1) % ny), (i, (j + 1) % ny)]


def boundary_fraction(labels: np.ndarray) -> float:
    """Fraction of unique nearest-neighbour bonds with unlike labels (periodic 2D)."""
    a = np.asarray(labels)
    if a.ndim != 2:
        raise ValueError("boundary_fraction expects a 2D array")
    right = np.roll(a, -1, axis=1)
    down = np.roll(a, -1, axis=0)
    return float((np.count_nonzero(a != right) + np.count_nonzero(a != down)) / (2 * a.size))
