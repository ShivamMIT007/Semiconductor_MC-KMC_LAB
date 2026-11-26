"""Classical Deal-Grove thermal oxidation baseline (continuum model)."""
from __future__ import annotations
import numpy as np
import pandas as pd
from semimc.common import make_result


def simulate(p: dict):
    A = float(p.get("A_nm", 120.0))
    B = float(p.get("B_nm2_s", 18.0))
    tau = float(p.get("initial_time_shift_s", 0.0))
    tmax = float(p.get("duration_s", 3600.0))
    n = int(p.get("points", 150))
    initial = float(p.get("initial_oxide_nm", 0.0))
    if A <= 0 or B <= 0 or tau < 0 or tmax <= 0 or n < 10 or initial < 0:
        raise ValueError("Deal-Grove parameters must be physically positive")
    t = np.linspace(0.0, tmax, n)
    # Solve x^2 + A*x = B*(t + tau_eff); choose tau_eff to honour x(0)=initial.
    tau_eff = tau + (initial * initial + A * initial) / B
    x = 0.5 * (-A + np.sqrt(A * A + 4.0 * B * (t + tau_eff)))
    df = pd.DataFrame({"time_s": t, "oxide_thickness_nm": x,
                       "instantaneous_growth_rate_nm_s": np.gradient(x, t)})
    before = np.zeros((40, 80), dtype=float)
    after = np.zeros_like(before)
    oxide_layers = max(1, min(39, int(round(x[-1] / max(x[-1], 1.0) * 20))))
    after[-oxide_layers:, :] = 1.0
    return make_result(
        name="thermal_oxidation", title="Silicon thermal oxidation · Deal–Grove model",
        method="Continuum analytical model · classical Deal–Grove baseline",
        description="Predicts oxide thickness from linear-parabolic growth kinetics; no stochastic events are used.",
        before=before, after=after, data=df,
        series={"Oxide thickness": {"x": t, "y": x, "xlabel": "Time (s)", "ylabel": "Oxide thickness (nm)"},
                 "Growth rate": {"x": t, "y": np.gradient(x, t), "xlabel": "Time (s)", "ylabel": "Growth rate (nm/s)"}},
        metrics={"final_oxide_thickness_nm": float(x[-1]),
                 "initial_oxide_thickness_nm": float(x[0]),
                 "linear_parabolic_A_nm": A, "parabolic_rate_B_nm2_s": B},
        visual_kind="cross_section", before_label="Bare silicon (schematic)", after_label="Oxidized silicon (schematic)",
        notes=["Deal–Grove is a continuum reference model and is deliberately not labelled MC or KMC.",
               "A and B are condition-dependent: temperature, ambient, pressure, orientation and stress matter.",
               "The cross-section image is schematic; thickness is represented quantitatively in the CSV and graph."],
        warnings=["A and B defaults are illustrative. Use measured or literature values for the exact oxidation ambient and temperature."],
        metadata={"A_nm": A, "B_nm2_s": B, "initial_time_shift_s": tau,
                  "duration_s": tmax, "initial_oxide_nm": initial},
    )
