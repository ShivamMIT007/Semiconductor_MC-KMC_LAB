from __future__ import annotations
import json
import streamlit as st
import pandas as pd
from semimc.catalog import CATALOG
from semimc.visuals import state_png, graphs_png, animation_gif, result_zip, save_result

st.set_page_config(page_title="Semiconductor MC/KMC Lab", page_icon="⚛️", layout="wide")

st.markdown("""
<style>
.block-container {padding-top: 1.4rem; padding-bottom: 2rem; max-width: 1500px;}
[data-testid="stMetric"] {background: rgba(100,116,139,.08); padding: 12px 14px; border-radius: 12px;}
.small-note {font-size: .88rem; color: #64748b;}
</style>
""", unsafe_allow_html=True)

st.title("⚛️ Semiconductor Monte Carlo / Kinetic Monte Carlo Lab")
st.caption("A modular Python simulation workbench for semiconductor processing and materials microstructure.")

with st.sidebar:
    st.header("Simulation selection")
    model_name = st.selectbox("Choose a process", list(CATALOG.keys()))
    st.divider()
    st.subheader("Model inputs")
    spec = CATALOG[model_name]
    params = {}
    with st.form("simulation_parameters"):
        for key, label, kind, default, minimum, maximum, step, help_text in spec["params"]:
            if kind == "select":
                params[key] = st.selectbox(label, ["B", "P", "As"], index=["B", "P", "As"].index(default), help=help_text)
            elif kind == "int":
                params[key] = st.number_input(label, min_value=int(minimum), max_value=int(maximum), value=int(default), step=int(step), help=help_text)
            else:
                params[key] = st.number_input(label, min_value=float(minimum), max_value=float(maximum), value=float(default), step=float(step), format="%.6g", help=help_text)
        run = st.form_submit_button("▶ Run simulation", type="primary", use_container_width=True)

    st.divider()
    st.markdown("**Model family guide**")
    st.markdown("- **MC:** random particle histories / configuration sampling\n- **KMC:** event rates + physical time\n- **Hybrid:** coupled transport and surface events\n- **Continuum:** deterministic reference physics")
    st.caption("Outputs are model demonstrations unless calibrated to material-specific data.")

if "last_result" not in st.session_state:
    st.session_state.last_result = None
if "last_model" not in st.session_state:
    st.session_state.last_model = None

if run:
    try:
        with st.spinner("Running stochastic simulation..."):
            result = spec["fn"](params)
        st.session_state.last_result = result
        st.session_state.last_model = model_name
        st.success(f"Completed: {result.title}")
    except Exception as exc:
        st.error(f"Simulation failed: {exc}")
        st.exception(exc)

result = st.session_state.last_result
if result is None:
    st.info("Choose a simulation in the left panel, adjust its inputs, and select **Run simulation**. The example-output folder contains reference runs generated with fixed seeds.")
    st.subheader("Available process models")
    cols = st.columns(3)
    descriptions = [
        ("Ion implantation", "MC stochastic range / straggle"),
        ("Dopant diffusion", "KMC on a 3D diamond-cubic network"),
        ("Thin-film deposition", "KMC solid-on-solid growth"),
        ("Grain growth", "Metropolis Potts MC"),
        ("Binary alloy ordering", "Conserved Kawasaki MC"),
        ("Phase separation", "Conserved Kawasaki MC"),
        ("Defect annealing", "Object KMC defect migration/recombination"),
        ("Plasma etching", "MC ion incidence / surface-yield prototype"),
        ("Thermal oxidation", "Deal–Grove continuum baseline"),
        ("Photolithography", "MC photon absorption"),
    ]
    for idx, (title, method) in enumerate(descriptions):
        with cols[idx % 3]:
            st.markdown(f"**{title}**")
            st.caption(method)
    st.stop()

st.header(result.title)
st.caption(result.method)
st.write(result.description)

metric_items = list(result.metrics.items())
if metric_items:
    display_items = metric_items[:4]
    metric_cols = st.columns(min(4, len(display_items)))
    for col, (key, value) in zip(metric_cols, display_items):
        if isinstance(value, float):
            text = f"{value:.4g}"
        else:
            text = str(value)
        col.metric(key.replace("_", " ").title(), text)

if result.warnings:
    for warning in result.warnings:
        st.warning(warning)

tabs = st.tabs(["Before / after", "Graphs", "Numerical dataset", "Method & validation", "Export"])
with tabs[0]:
    c1, c2 = st.columns(2)
    c1.image(state_png(result, "before"), caption=result.before_label, use_container_width=True)
    c2.image(state_png(result, "after"), caption=result.after_label, use_container_width=True)
    gif = animation_gif(result)
    if gif:
        st.markdown("**Evolution animation**")
        st.image(gif, caption="Saved simulation snapshots", use_container_width=True)
        st.download_button("Download evolution GIF", gif, file_name=f"{result.name}_evolution.gif", mime="image/gif")
    if result.snapshots:
        st.caption(f"Stored {len(result.snapshots)} snapshots. The animation shows selected saved states, not every event.")

with tabs[1]:
    st.image(graphs_png(result), caption="Simulation time-series / output plots", use_container_width=True)
    if result.series:
        for title, series in result.series.items():
            st.markdown(f"**{title}**")
            chart_df = pd.DataFrame({series.get("xlabel", "x"): series["x"], series.get("ylabel", "y"): series["y"]})
            st.line_chart(chart_df.set_index(chart_df.columns[0]))

with tabs[2]:
    st.markdown(f"**Rows:** {len(result.data):,} · **Columns:** {len(result.data.columns)}")
    st.dataframe(result.data, use_container_width=True, height=420)
    st.download_button("Download numerical dataset (CSV)", result.data.to_csv(index=False).encode("utf-8"),
                       file_name=f"{result.name}_dataset.csv", mime="text/csv")

with tabs[3]:
    st.subheader("Model notes and assumptions")
    for item in result.notes:
        st.markdown(f"- {item}")
    st.subheader("Numerical output metrics")
    st.json(result.metrics)
    st.subheader("Simulation metadata")
    st.json(result.metadata or {})
    st.info("Research validation requires algorithm verification, conservation checks, convergence studies, material-specific parameter provenance, and comparison with analytical or experimental reference data. A visually plausible image is not validation.")

with tabs[4]:
    st.markdown("Export this run as a self-contained research record.")
    try:
        bundle = result_zip(result)
        st.download_button("Download complete run bundle (ZIP)", bundle,
                           file_name=f"{result.name}_run_bundle.zip", mime="application/zip", type="primary")
    except Exception as exc:
        st.error(f"Could not create run bundle: {exc}")
    if st.button("Save this run to local results/ folder"):
        try:
            folder = save_result(result)
            st.success(f"Saved files to: {folder}")
        except Exception as exc:
            st.error(f"Could not save run: {exc}")

st.divider()
st.caption("Semiconductor MC/KMC Lab · Research prototype. Confirm model assumptions and parameter sources before using results in publications or process decisions.")
