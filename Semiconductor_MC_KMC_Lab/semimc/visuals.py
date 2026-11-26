from __future__ import annotations
from io import BytesIO
import json
import zipfile
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import colormaps
from PIL import Image
from semimc.common import SimulationResult


def _state_image(arr: np.ndarray, kind: str, label: str) -> bytes:
    a = np.asarray(arr)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    if a.ndim == 1 or (a.ndim == 2 and min(a.shape) == 1):
        y = a.reshape(-1)
        ax.plot(np.arange(y.size), y, linewidth=1.8)
        ax.set_xlabel("Lateral position / index")
        ax.set_ylabel("Value")
        ax.grid(alpha=0.25)
    else:
        if kind == "labels":
            cmap = colormaps["tab20"]
            im = ax.imshow(a, cmap=cmap, interpolation="nearest", origin="lower")
            fig.colorbar(im, ax=ax, label="Grain label")
        elif kind == "spins":
            im = ax.imshow(a, cmap="coolwarm", vmin=-1, vmax=1, interpolation="nearest", origin="lower")
            fig.colorbar(im, ax=ax, ticks=[-1, 1], label="Species / spin")
        elif kind == "defects":
            cmap = ListedColormap(["#F4F4F4", "#E45756", "#4C78A8"])
            norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
            im = ax.imshow(a, cmap=cmap, norm=norm, interpolation="nearest", origin="lower")
            cb = fig.colorbar(im, ax=ax, ticks=[0, 1, 2])
            cb.ax.set_yticklabels(["Empty", "Vacancy", "Interstitial"])
        elif kind in {"height", "heatmap", "implant"}:
            cmap = "viridis" if kind == "height" else "magma"
            im = ax.imshow(a, cmap=cmap, interpolation="nearest", origin="upper" if kind == "implant" else "lower", aspect="auto")
            fig.colorbar(im, ax=ax, label="Height / intensity / count")
        else:
            im = ax.imshow(a, cmap="viridis", interpolation="nearest", origin="lower", aspect="auto")
            fig.colorbar(im, ax=ax, label="State value")
        ax.set_xlabel("Lateral lattice index")
        ax.set_ylabel("Depth / lattice index")
    ax.set_title(label)
    fig.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def state_png(result: SimulationResult, which: str = "after") -> bytes:
    arr = result.before if which == "before" else result.after
    label = result.before_label if which == "before" else result.after_label
    return _state_image(arr, result.visual_kind, label)


def graphs_png(result: SimulationResult) -> bytes:
    items = list(result.series.items())
    if not items:
        fig, ax = plt.subplots(figsize=(7, 3), dpi=140)
        ax.text(0.5, 0.5, "No time-series data for this model", ha="center", va="center")
        ax.axis("off")
    else:
        fig, axes = plt.subplots(len(items), 1, figsize=(8, max(3, 2.8 * len(items))), dpi=140, squeeze=False)
        for ax, (title, s) in zip(axes[:, 0], items):
            ax.plot(np.asarray(s["x"]), np.asarray(s["y"]), linewidth=1.7)
            ax.set_title(title)
            ax.set_xlabel(s.get("xlabel", "x"))
            ax.set_ylabel(s.get("ylabel", "y"))
            ax.grid(alpha=0.25)
        fig.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def animation_gif(result: SimulationResult, duration_ms: int = 350) -> bytes | None:
    snapshots = result.snapshots or []
    if len(snapshots) < 2:
        return None
    frames = []
    kind = result.visual_kind
    stack = np.stack([np.asarray(s) for s in snapshots])
    vmin, vmax = float(np.min(stack)), float(np.max(stack))
    if kind == "labels": cmap_name = "tab20"
    elif kind == "spins": cmap_name = "coolwarm"
    elif kind == "height": cmap_name = "viridis"
    elif kind == "defects": cmap_name = "viridis"
    else: cmap_name = "magma"
    cmap = colormaps[cmap_name]
    for arr in snapshots:
        a = np.asarray(arr, dtype=float)
        if vmax == vmin:
            norm = np.zeros_like(a)
        else:
            norm = (a - vmin) / (vmax - vmin)
        rgba = (cmap(norm)[:, :, :3] * 255).astype(np.uint8)
        img = Image.fromarray(rgba, mode="RGB").resize((480, 360), Image.Resampling.NEAREST)
        frames.append(img)
    buf = BytesIO()
    frames[0].save(buf, format="GIF", save_all=True, append_images=frames[1:],
                   duration=duration_ms, loop=0, optimize=False)
    return buf.getvalue()


def result_zip(result: SimulationResult) -> bytes:
    buf = BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("summary.json", json.dumps(result.summary(), indent=2, ensure_ascii=False))
        zf.writestr("dataset.csv", result.data.to_csv(index=False))
        zf.writestr("before.png", state_png(result, "before"))
        zf.writestr("after.png", state_png(result, "after"))
        zf.writestr("graphs.png", graphs_png(result))
        gif = animation_gif(result)
        if gif:
            zf.writestr("evolution.gif", gif)
        if result.metadata:
            safe_meta = {k: v for k, v in result.metadata.items() if isinstance(v, (str, int, float, bool, list, dict, type(None)))}
            zf.writestr("input_parameters.json", json.dumps(safe_meta, indent=2, ensure_ascii=False, default=str))
    return buf.getvalue()


def save_result(result: SimulationResult, base_dir: str = "results") -> Path:
    root = Path(base_dir)
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = root / f"{result.name}_{stamp}"
    out.mkdir(parents=True, exist_ok=False)
    (out / "summary.json").write_text(json.dumps(result.summary(), indent=2, ensure_ascii=False), encoding="utf-8")
    result.data.to_csv(out / "dataset.csv", index=False)
    (out / "before.png").write_bytes(state_png(result, "before"))
    (out / "after.png").write_bytes(state_png(result, "after"))
    (out / "graphs.png").write_bytes(graphs_png(result))
    gif = animation_gif(result)
    if gif:
        (out / "evolution.gif").write_bytes(gif)
    return out
