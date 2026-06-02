"""Streamlit web app for a PubMed-backed network meta-analysis assistant."""

from __future__ import annotations

import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from pubmed_agent import (
    articles_to_dataframe,
    build_evidence_network,
    build_pubmed_query,
    generate_protocol_notes,
    infer_treatment_mentions,
    normalize_terms,
    pubmed_search_url,
    retrieve_references,
)

st.set_page_config(page_title="PubMed NMA Agent", page_icon="🧬", layout="wide")

st.title("🧬 PubMed Network Meta-analysis Agent")
st.caption(
    "Build a transparent PICO search, retrieve PubMed-indexed references, and draft the evidence network for a network meta-analysis."
)

with st.sidebar:
    st.header("Search configuration")
    condition = st.text_input("Condition / population", placeholder="e.g. rheumatoid arthritis")
    interventions = st.text_area("Interventions (comma, semicolon, or newline separated)", placeholder="adalimumab\netanercept\ninfliximab")
    comparators = st.text_area("Comparators / placebo / standard care", placeholder="placebo\nmethotrexate")
    outcome = st.text_input("Outcome", placeholder="ACR20 response")
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
st.write("After screening full texts, enter direct comparisons found in eligible studies. Each row represents one study comparison.")
default_edges = pd.DataFrame(
    {
        "study_id": [""],
        "treatment_a": [""],
        "treatment_b": [""],
        "pmid": [""],
    }
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
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(width=1.5), hoverinfo="none"))
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=list(network.nodes()),
            textposition="top center",
            marker=dict(size=22, color="#2E86AB"),
        )
    )
    fig.update_layout(showlegend=False, margin=dict(l=20, r=20, t=20, b=20))
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Add at least one valid treatment comparison to visualize the network.")

st.subheader("5. NMA protocol checklist")
for note in generate_protocol_notes(condition, treatment_terms, outcome):
    st.markdown(f"- {note}")

st.warning(
    "This app supports literature discovery and planning; it is not a substitute for duplicate screening, full-text data extraction, "
    "clinical/statistical review, or a validated NMA model implementation."
)
