"""Streamlit web app for a PubMed-backed network meta-analysis assistant."""

from __future__ import annotations

import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from pubmed_agent import (
    DEFAULT_TOPIC,
    articles_to_dataframe,
    build_evidence_network,
    build_pubmed_query,
    dataframe_to_tiff,
    demo_extracted_comparisons,
    generate_protocol_notes,
    infer_treatment_mentions,
    normalize_terms,
    plotly_figure_to_tiff,
    pubmed_search_url,
    retrieve_references,
    summarize_network,
)

st.set_page_config(page_title="PubMed NMA Agent", page_icon="🧬", layout="wide")

st.title("🧬 PubMed Network Meta-analysis Agent")
st.caption(
    "Build a transparent PICO search, retrieve PubMed-indexed references, and draft the evidence network for a network meta-analysis."
)

with st.sidebar:
    st.header("Search configuration")
    st.write("Default topic: biologic/DMARD therapies for rheumatoid arthritis and ACR20 response.")
    condition = st.text_input("Condition / population", value=DEFAULT_TOPIC["condition"])
    interventions = st.text_area("Interventions (comma, semicolon, or newline separated)", value=DEFAULT_TOPIC["interventions"])
    comparators = st.text_area("Comparators / placebo / standard care", value=DEFAULT_TOPIC["comparators"])
    outcome = st.text_input("Outcome", value=DEFAULT_TOPIC["outcome"])
    max_results = st.slider("Maximum PubMed references", min_value=5, max_value=100, value=25, step=5)
    include_rct_filter = st.checkbox("Prioritize randomized trials", value=True)
    include_nma_filter = st.checkbox("Search existing NMAs / indirect comparisons", value=False)
    include_humans_english = st.checkbox("Humans + English filter", value=True)
    email = st.text_input("NCBI email (recommended)", placeholder="you@example.org")
    run_search = st.button("Search PubMed", type="primary")

treatment_terms = normalize_terms(interventions) + normalize_terms(comparators)

query = build_pubmed_query(
    condition=condition,
    interventions=interventions,
    comparators=comparators,
    outcome=outcome,
    include_rct_filter=include_rct_filter,
    include_nma_filter=include_nma_filter,
    include_humans_english=include_humans_english,
)

st.subheader("1. Reproducible PubMed query")
st.code(query or "Enter a condition and interventions to build a query.", language="text")
if query:
    st.link_button("Open this query in PubMed", pubmed_search_url(query))

st.subheader("2. PubMed references")
if run_search:
    if not query:
        st.warning("Please provide enough PICO terms to build a PubMed query.")
    else:
        with st.spinner("Searching PubMed via NCBI E-utilities..."):
            try:
                st.session_state["articles"] = retrieve_references(query, max_results=max_results, email=email)
            except Exception as exc:  # Streamlit should show network/API errors to the user.
                st.error(f"PubMed retrieval failed: {exc}")
                st.session_state["articles"] = []

articles = st.session_state.get("articles", [])
if articles:
    refs = articles_to_dataframe(articles)
    st.success(f"Retrieved {len(refs)} PubMed-indexed references.")
    st.dataframe(
        refs[["pmid", "year", "title", "journal", "authors", "publication_types", "url"]],
        hide_index=True,
        use_container_width=True,
    )
    st.download_button(
        "Download references as CSV",
        refs.to_csv(index=False).encode("utf-8"),
        file_name="pubmed_nma_references.csv",
        mime="text/csv",
    )

    st.subheader("3. AI-assisted screening cues")
    mentions = infer_treatment_mentions(articles, treatment_terms)
    st.dataframe(mentions, hide_index=True, use_container_width=True)
    st.info("Treatment mentions are keyword-based screening cues. Confirm study arms and extracted outcomes from full-text reports before analysis.")
else:
    st.info("Run a PubMed search to populate the reference table and screening cues.")

st.subheader("4. Draft treatment network")
st.write(
    "After screening full texts, enter direct comparisons found in eligible studies. "
    "For the NMA summary, enter log effect estimates (for example log odds ratios) and standard errors; higher values favour treatment A."
)
use_demo = st.toggle("Load synthetic demonstration contrasts for the default topic", value=True)
default_edges = demo_extracted_comparisons() if use_demo else pd.DataFrame(
    {"study_id": [""], "treatment_a": [""], "treatment_b": [""], "effect": [""], "se": [""], "pmid": [""]}
)
edges = st.data_editor(default_edges, num_rows="dynamic", use_container_width=True, key="edge_editor")
network = build_evidence_network(edges)
if network.number_of_edges() > 0:
    positions = nx.spring_layout(network, seed=42)
    edge_x: list[float | None] = []
    edge_y: list[float | None] = []
    for start, end in network.edges():
        x0, y0 = positions[start]
        x1, y1 = positions[end]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    node_x = [positions[node][0] for node in network.nodes()]
    node_y = [positions[node][1] for node in network.nodes()]
    fig = go.Figure()
    edge_weights = [network[u][v].get("weight", 1) for u, v in network.edges()]
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(width=1.5), hoverinfo="none"))
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=list(network.nodes()),
            textposition="top center",
            marker=dict(size=28, color="#2E86AB"),
        )
    )
    fig.update_layout(title=f"Treatment evidence network ({sum(edge_weights)} direct comparison row[s])", showlegend=False, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(fig, use_container_width=True)
    try:
        st.download_button(
            "Download network plot as TIFF",
            plotly_figure_to_tiff(fig),
            file_name="nma_network.tif",
            mime="image/tiff",
        )
    except Exception:
        st.caption("TIFF export for Plotly figures requires the optional Kaleido renderer in the runtime.")
else:
    st.info("Add at least one valid treatment comparison to visualize the network.")

st.subheader("5. Exploratory NMA summary")
candidate_refs = sorted(network.nodes()) if network.number_of_nodes() else normalize_terms(comparators) + normalize_terms(interventions)
reference = st.selectbox("Reference treatment", candidate_refs or ["placebo"], index=0)
summary = summarize_network(edges, reference=reference)
st.markdown(summary["narrative"])
if not summary["relative"].empty:
    relative = summary["relative"].copy()
    relative["effect (95% CI)"] = relative.apply(lambda row: f"{row['effect']:.2f} ({row['ci_low']:.2f}, {row['ci_high']:.2f})", axis=1)
    relative["ratio"] = relative["ratio"].map(lambda value: f"{value:.2f}")
    relative["rank_score"] = relative["rank_score"].map(lambda value: f"{value:.2f}")
    st.dataframe(relative[["treatment", "effect (95% CI)", "ratio", "rank_score", "paths"]], hide_index=True, use_container_width=True)

    forest = go.Figure()
    forest.add_trace(
        go.Scatter(
            x=summary["relative"]["effect"],
            y=summary["relative"]["treatment"],
            error_x=dict(
                type="data",
                symmetric=False,
                array=summary["relative"]["ci_high"] - summary["relative"]["effect"],
                arrayminus=summary["relative"]["effect"] - summary["relative"]["ci_low"],
            ),
            mode="markers",
            marker=dict(size=10, color="#D1495B"),
        )
    )
    forest.add_vline(x=0, line_dash="dash", line_color="gray")
    forest.update_layout(title=f"Relative effects versus {summary['reference']}", xaxis_title="Log relative effect", yaxis_title="", margin=dict(l=20, r=20, t=50, b=40))
    st.plotly_chart(forest, use_container_width=True)
    try:
        st.download_button("Download forest plot as TIFF", plotly_figure_to_tiff(forest), file_name="nma_forest.tif", mime="image/tiff")
    except Exception:
        st.caption("TIFF export for Plotly figures requires the optional Kaleido renderer in the runtime.")

    st.write("League table (row treatment versus column treatment; ratio scale)")
    st.dataframe(summary["league"], hide_index=True, use_container_width=True)
    st.download_button(
        "Download league table as TIFF",
        dataframe_to_tiff(summary["league"], title="NMA league table"),
        file_name="nma_league_table.tif",
        mime="image/tiff",
    )
else:
    st.info("Add effect estimates and standard errors to generate the exploratory summary.")

st.subheader("6. NMA protocol checklist")
for note in generate_protocol_notes(condition, treatment_terms, outcome):
    st.markdown(f"- {note}")

st.warning(
    "This app supports literature discovery and planning; it is not a substitute for duplicate screening, full-text data extraction, "
    "clinical/statistical review, or a validated NMA model implementation."
)
