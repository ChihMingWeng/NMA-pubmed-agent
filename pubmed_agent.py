"""Utilities for a PubMed-backed network meta-analysis assistant.

The module intentionally avoids doing statistical NMA itself.  It helps users
build reproducible PubMed searches, retrieve PubMed-indexed references, and
summarize the evidence network that should be exported to a specialist package
such as R netmeta/gemtc, BUGS/JAGS, or Python workflows.
"""

from __future__ import annotations

import html
import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from io import BytesIO
from math import erf, exp, sqrt
from typing import Any, Iterable
from urllib.parse import urlencode


NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED_ARTICLE_URL = "https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
DEFAULT_TOOL = "nma-pubmed-agent"
DEFAULT_TOPIC = {
    "condition": "severe asthma",
    "interventions": "dupilumab\nomalizumab\nmepolizumab\nbenralizumab\ntezepelumab\nreslizumab",
    "comparators": "placebo\nstandard care",
    "outcome": "exacerbation",
}
DISEASE_TEMPLATES = {
    "Rheumatoid Arthritis": {
        "condition": "rheumatoid arthritis",
        "interventions": "adalimumab\netanercept\ninfliximab\nmethotrexate",
        "comparators": "placebo\nstandard care",
        "outcome": "ACR20 response",
    },
    "Severe Asthma Biologics": DEFAULT_TOPIC,
    "Asthma and COVID-19": {
        "condition": "asthma; COVID-19",
        "interventions": "omalizumab\nmepolizumab\nbenralizumab\ndupilumab\ntezepelumab\nreslizumab\nICS\nICS/LABA\noral corticosteroids",
        "comparators": "no biologic treatment\nstandard care\nplacebo",
        "outcome": "hospitalization\nICU admission\nmechanical ventilation\nmortality",
    },
}
STUDY_DESIGNS = ("Randomized Controlled Trials", "Observational Studies", "Both")

RCT_FILTER = '((randomized controlled trial[Publication Type]) OR randomized[Title/Abstract] OR randomised[Title/Abstract] OR trial[Title/Abstract])'
OBSERVATIONAL_FILTER = (
    "((observational study[Publication Type]) OR cohort[Title/Abstract] OR case-control[Title/Abstract] "
    "OR cross-sectional[Title/Abstract] OR registry[Title/Abstract] OR real-world[Title/Abstract])"
)
NMA_FILTER = '((network meta-analysis[Title/Abstract]) OR mixed treatment comparison[Title/Abstract] OR indirect treatment comparison[Title/Abstract])'
HUMANS_ENGLISH_FILTER = '(humans[MeSH Terms] AND english[Language])'

VALIDATED_CONDITION_MESH = {
    "asthma": "Asthma",
    "covid-19": "COVID-19",
    "rheumatoid arthritis": "Arthritis, Rheumatoid",
}
BIOLOGIC_TERMS = {
    "dupilumab",
    "omalizumab",
    "mepolizumab",
    "benralizumab",
    "tezepelumab",
    "reslizumab",
}
PLACEBO_TERMS = {"placebo", "placebo control", "placebo controlled", "placebo-controlled"}


@dataclass(frozen=True)
class PubMedArticle:
    """Structured data for a PubMed-indexed article."""

    pmid: str
    title: str
    journal: str
    year: str
    authors: str
    abstract: str
    url: str
    publication_types: str


def normalize_terms(raw_terms: str | Iterable[str]) -> list[str]:
    """Return non-empty search terms with repeated whitespace collapsed."""

    if isinstance(raw_terms, str):
        parts = re.split(r"[,;\n]", raw_terms)
    else:
        parts = list(raw_terms)
    normalized = []
    for part in parts:
        if not part or not part.strip():
            continue
        term = re.sub(r"\s+", " ", part).strip()
        if term.lower() in {"placebo/treatment", "placebo treatment"}:
            term = "placebo"
        normalized.append(term)
    return normalized


def _quote_term(term: str) -> str:
    return term.replace('"', "")


def _fielded_or(terms: Iterable[str], fields: tuple[str, ...]) -> str:
    clauses: list[str] = []
    for term in normalize_terms(terms):
        escaped = _quote_term(term)
        field_clauses = [f'"{escaped}"[{field}]' for field in fields]
        clauses.append("(" + " OR ".join(field_clauses) + ")")
    return "(" + " OR ".join(clauses) + ")" if clauses else ""


def _condition_clause(condition: str) -> str:
    """Build a condition clause with curated MeSH only for validated headings."""

    terms = normalize_terms(condition)
    clauses = []
    for term in terms:
        escaped = _quote_term(term)
        title_abstract = f'"{escaped}"[Title/Abstract]'
        mesh = VALIDATED_CONDITION_MESH.get(term.lower())
        if mesh and mesh.lower() == term.lower():
            clauses.append(f'("{_quote_term(mesh)}"[MeSH Terms] OR {title_abstract})')
        elif mesh:
            clauses.append(f"({title_abstract} OR \"{_quote_term(mesh)}\"[MeSH Terms])")
        else:
            clauses.append(f"({title_abstract})")
    return "(" + " AND ".join(clauses) + ")" if clauses else ""


def _study_design_clause(study_design: str) -> str:
    """Return a PubMed study-design clause for the selected design."""

    if study_design == "Randomized Controlled Trials":
        return RCT_FILTER
    if study_design == "Observational Studies":
        return OBSERVATIONAL_FILTER
    if study_design == "Both":
        return f"({RCT_FILTER} OR {OBSERVATIONAL_FILTER})"
    return ""


def _treatment_clause(interventions: str | Iterable[str], comparators: str | Iterable[str]) -> str:
    """Build intervention/comparator clauses without assuming drug MeSH headings exist."""

    clauses = []
    for term in normalize_terms(interventions):
        escaped = _quote_term(term)
        if term.lower() in BIOLOGIC_TERMS:
            clauses.append(f'("{escaped}"[Title/Abstract] OR "{escaped}"[Substance Name])')
        else:
            clauses.append(f'("{escaped}"[Title/Abstract])')
    for term in normalize_terms(comparators):
        escaped = _quote_term(term)
        if term.lower() in PLACEBO_TERMS:
            clauses.append('("placebo"[Title/Abstract] OR "Placebos"[MeSH Terms])')
        else:
            clauses.append(f'("{escaped}"[Title/Abstract])')
    return "(" + " OR ".join(clauses) + ")" if clauses else ""


def build_pubmed_query(
    condition: str,
    interventions: str | Iterable[str],
    comparators: str | Iterable[str] = "",
    outcome: str = "",
    include_rct_filter: bool = True,
    include_nma_filter: bool = False,
    include_humans_english: bool = True,
    study_design: str | None = None,
) -> str:
    """Build a transparent PubMed query from a PICO-style question."""

    clauses: list[str] = []
    condition_clause = _condition_clause(condition) if condition else ""
    treatment_clause = _treatment_clause(interventions, comparators)
    outcome_clause = _fielded_or(outcome, ("Title/Abstract",)) if outcome else ""

    for clause in (condition_clause, treatment_clause, outcome_clause):
        if clause:
            clauses.append(clause)
    design_clause = _study_design_clause(study_design) if study_design else (RCT_FILTER if include_rct_filter else "")
    if design_clause:
        clauses.append(design_clause)
    if include_nma_filter:
        clauses.append(NMA_FILTER)
    if include_humans_english:
        clauses.append(HUMANS_ENGLISH_FILTER)
    return " AND ".join(clauses)


def pubmed_search_url(query: str) -> str:
    """Return the PubMed web URL for a query."""

    return f"https://pubmed.ncbi.nlm.nih.gov/?{urlencode({'term': query})}"


def search_pubmed(query: str, max_results: int = 25, email: str = "") -> list[str]:
    """Search PubMed and return PMIDs ordered by PubMed relevance/date settings."""

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max(1, min(max_results, 200)),
        "tool": DEFAULT_TOOL,
        "sort": "relevance",
    }
    if email:
        params["email"] = email
    import requests

    response = requests.get(f"{NCBI_EUTILS}/esearch.fcgi", params=params, timeout=20)
    response.raise_for_status()
    payload = response.json()
    return payload.get("esearchresult", {}).get("idlist", [])


def fetch_pubmed_details(pmids: Iterable[str], email: str = "") -> list[PubMedArticle]:
    """Fetch PubMed article details for PMIDs using NCBI EFetch."""

    ids = [pmid for pmid in pmids if pmid]
    if not ids:
        return []
    params = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "xml",
        "tool": DEFAULT_TOOL,
    }
    if email:
        params["email"] = email
    import requests

    response = requests.get(f"{NCBI_EUTILS}/efetch.fcgi", params=params, timeout=30)
    response.raise_for_status()
    root = ET.fromstring(response.content)
    return [_parse_article(article) for article in root.findall(".//PubmedArticle")]


def _node_text(node: ET.Element | None, default: str = "") -> str:
    if node is None:
        return default
    return html.unescape(" ".join("".join(node.itertext()).split()))


def _parse_article(article: ET.Element) -> PubMedArticle:
    medline = article.find("MedlineCitation")
    article_node = medline.find("Article") if medline is not None else None
    pmid = _node_text(medline.find("PMID") if medline is not None else None)
    title = _node_text(article_node.find("ArticleTitle") if article_node is not None else None, "Untitled")
    journal = _node_text(article_node.find("Journal/Title") if article_node is not None else None, "Unknown journal")
    year = _node_text(article_node.find("Journal/JournalIssue/PubDate/Year") if article_node is not None else None)
    if not year:
        medline_date = _node_text(article_node.find("Journal/JournalIssue/PubDate/MedlineDate") if article_node is not None else None)
        year_match = re.search(r"\d{4}", medline_date)
        year = year_match.group(0) if year_match else "Unknown year"
    author_nodes = article_node.findall("AuthorList/Author") if article_node is not None else []
    author_names = []
    for author in author_nodes[:6]:
        last = _node_text(author.find("LastName"))
        initials = _node_text(author.find("Initials"))
        collective = _node_text(author.find("CollectiveName"))
        if collective:
            author_names.append(collective)
        elif last:
            author_names.append(f"{last} {initials}".strip())
    authors = ", ".join(author_names)
    if len(author_nodes) > 6:
        authors += ", et al."
    abstract_parts = article_node.findall("Abstract/AbstractText") if article_node is not None else []
    abstract = " ".join(_node_text(part) for part in abstract_parts) or "No abstract available from PubMed."
    pub_types = article_node.findall("PublicationTypeList/PublicationType") if article_node is not None else []
    publication_types = "; ".join(_node_text(pub_type) for pub_type in pub_types)
    return PubMedArticle(
        pmid=pmid,
        title=title,
        journal=journal,
        year=year,
        authors=authors or "Authors not listed",
        abstract=abstract,
        url=PUBMED_ARTICLE_URL.format(pmid=pmid),
        publication_types=publication_types,
    )


def retrieve_references(query: str, max_results: int = 25, email: str = "") -> list[PubMedArticle]:
    """Search PubMed and fetch article metadata, respecting NCBI rate limits."""

    pmids = search_pubmed(query, max_results=max_results, email=email)
    time.sleep(0.34)
    return fetch_pubmed_details(pmids, email=email)


def articles_to_dataframe(articles: Iterable[PubMedArticle]) -> Any:
    """Convert articles to a table suitable for Streamlit display/export."""

    import pandas as pd

    return pd.DataFrame([article.__dict__ for article in articles])


def nma_dataset_from_edges(edge_rows: Any) -> Any:
    """Return an export-ready NMA dataset table from edited comparison rows."""

    import pandas as pd

    columns = {
        "pmid": "PMID",
        "treatment_a": "treatment",
        "treatment_b": "comparator",
        "effect": "effect estimate",
        "se": "standard error",
    }
    if edge_rows.empty:
        return pd.DataFrame(columns=list(columns.values()))
    dataset = edge_rows.copy()
    for source in columns:
        if source not in dataset.columns:
            dataset[source] = ""
    dataset = dataset[list(columns)].rename(columns=columns)
    return dataset


def infer_treatment_mentions(articles: Iterable[PubMedArticle], treatments: Iterable[str]) -> Any:
    """Flag article-level treatment mentions in title/abstract text.

    This is a screening aid only; users must verify study arms manually before
    performing an NMA.
    """

    import pandas as pd

    rows = []
    normalized_treatments = normalize_terms(treatments)
    for article in articles:
        haystack = f"{article.title} {article.abstract}".lower()
        mentions = [term for term in normalized_treatments if term.lower() in haystack]
        rows.append(
            {
                "pmid": article.pmid,
                "title": article.title,
                "mentioned_treatments": ", ".join(mentions),
                "mention_count": len(mentions),
                "url": article.url,
            }
        )
    return pd.DataFrame(rows).sort_values(["mention_count", "pmid"], ascending=[False, True]) if rows else pd.DataFrame()


def build_evidence_network(edge_rows: Any) -> Any:
    """Build a treatment network from rows with treatment_a/treatment_b columns."""

    import networkx as nx

    graph = nx.Graph()
    required = {"treatment_a", "treatment_b"}
    if edge_rows.empty or not required.issubset(edge_rows.columns):
        return graph
    for _, row in edge_rows.dropna(subset=["treatment_a", "treatment_b"]).iterrows():
        a = str(row["treatment_a"]).strip()
        b = str(row["treatment_b"]).strip()
        if not a or not b or a == b:
            continue
        weight = graph[a][b]["weight"] + 1 if graph.has_edge(a, b) else 1
        graph.add_edge(a, b, weight=weight)
    return graph


def demo_extracted_comparisons() -> Any:
    """Return example extracted contrast data for the built-in severe asthma topic.

    The rows are synthetic demonstration effect estimates intended to make the
    web app runnable before users finish full-text data extraction. Users should
    replace them with verified study-level contrasts before drawing conclusions.
    """

    import pandas as pd

    return pd.DataFrame(
        [
            {"study_id": "DEMO-ASTHMA-01", "treatment_a": "dupilumab", "treatment_b": "placebo", "effect": 0.63, "se": 0.18, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-02", "treatment_a": "omalizumab", "treatment_b": "placebo", "effect": 0.41, "se": 0.20, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-03", "treatment_a": "mepolizumab", "treatment_b": "placebo", "effect": 0.55, "se": 0.19, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-04", "treatment_a": "benralizumab", "treatment_b": "placebo", "effect": 0.50, "se": 0.21, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-05", "treatment_a": "tezepelumab", "treatment_b": "placebo", "effect": 0.69, "se": 0.17, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-06", "treatment_a": "reslizumab", "treatment_b": "placebo", "effect": 0.38, "se": 0.23, "pmid": ""},
            {"study_id": "DEMO-ASTHMA-07", "treatment_a": "dupilumab", "treatment_b": "standard care", "effect": 0.44, "se": 0.24, "pmid": ""},
        ]
    )


def summarize_network(edge_rows: Any, reference: str = "") -> dict[str, Any]:
    """Create simple NMA planning summaries from extracted log-effect contrasts.

    Effects are interpreted as log odds/risk/hazard ratios where higher values
    favour ``treatment_a`` over ``treatment_b``. The function pools duplicate
    direct contrasts by inverse variance and estimates each treatment versus the
    selected reference by the precision-weighted average of all simple paths.
    It is a transparent planning approximation, not a replacement for a
    validated NMA package.
    """

    import pandas as pd

    required = {"treatment_a", "treatment_b", "effect", "se"}
    if edge_rows.empty or not required.issubset(edge_rows.columns):
        return {"reference": reference, "direct": pd.DataFrame(), "relative": pd.DataFrame(), "league": pd.DataFrame(), "narrative": "Enter contrasts with effect and SE values to generate an NMA summary."}

    clean = edge_rows.copy()
    clean["treatment_a"] = clean["treatment_a"].astype(str).str.strip()
    clean["treatment_b"] = clean["treatment_b"].astype(str).str.strip()
    clean["effect"] = pd.to_numeric(clean["effect"], errors="coerce")
    clean["se"] = pd.to_numeric(clean["se"], errors="coerce")
    clean = clean[(clean["treatment_a"] != "") & (clean["treatment_b"] != "") & (clean["treatment_a"] != clean["treatment_b"]) & clean["effect"].notna() & clean["se"].gt(0)]
    if clean.empty:
        return {"reference": reference, "direct": pd.DataFrame(), "relative": pd.DataFrame(), "league": pd.DataFrame(), "narrative": "No valid contrasts are available yet."}

    direct_rows = []
    graph: dict[str, dict[str, tuple[float, float]]] = {}
    for (a, b), group in clean.groupby(["treatment_a", "treatment_b"], sort=True):
        weights = 1 / (group["se"] ** 2)
        pooled = float((group["effect"] * weights).sum() / weights.sum())
        se = float(sqrt(1 / weights.sum()))
        direct_rows.append({"treatment_a": a, "treatment_b": b, "effect": pooled, "se": se, "studies": int(len(group))})
        graph.setdefault(a, {})[b] = (pooled, se)
        graph.setdefault(b, {})[a] = (-pooled, se)
    direct = pd.DataFrame(direct_rows)
    treatments = sorted(set(clean["treatment_a"]).union(clean["treatment_b"]))
    reference = reference if reference in treatments else treatments[0]

    def path_estimates(start: str, end: str) -> list[tuple[float, float]]:
        estimates: list[tuple[float, float]] = []
        stack = [(start, [], 0.0, 0.0)]
        while stack:
            node, visited, estimate, variance = stack.pop()
            if node == end and visited:
                estimates.append((estimate, sqrt(variance)))
                continue
            if node in visited:
                continue
            for nxt, (edge_effect, edge_se) in graph.get(node, {}).items():
                if nxt not in visited:
                    stack.append((nxt, visited + [node], estimate + edge_effect, variance + edge_se**2))
        return estimates

    rel_rows = []
    for treatment in treatments:
        if treatment == reference:
            rel_rows.append({"treatment": treatment, "effect": 0.0, "se": 0.0, "ci_low": 0.0, "ci_high": 0.0, "ratio": 1.0, "rank_score": 0.5, "paths": 0})
            continue
        estimates = path_estimates(treatment, reference)
        weights = [1 / (se**2) for _, se in estimates if se > 0]
        if not estimates or not weights:
            continue
        effect = sum(est * wt for (est, se), wt in zip(estimates, weights) if se > 0) / sum(weights)
        se = sqrt(1 / sum(weights))
        z = effect / se if se else 0.0
        rank_score = 0.5 * (1 + erf(z / sqrt(2)))
        rel_rows.append({"treatment": treatment, "effect": effect, "se": se, "ci_low": effect - 1.96 * se, "ci_high": effect + 1.96 * se, "ratio": exp(effect), "rank_score": rank_score, "paths": len(estimates)})
    relative = pd.DataFrame(rel_rows).sort_values("effect", ascending=False)

    league = pd.DataFrame(index=treatments, columns=treatments, dtype=object)
    rel_map = {row["treatment"]: row for _, row in relative.iterrows()}
    for row in treatments:
        for col in treatments:
            if row == col:
                league.loc[row, col] = "—"
            elif row in rel_map and col in rel_map:
                diff = rel_map[row]["effect"] - rel_map[col]["effect"]
                ratio = exp(diff)
                league.loc[row, col] = f"{ratio:.2f}"
    best = relative.iloc[0]["treatment"] if not relative.empty else "No treatment"
    narrative = (
        f"Using {reference} as the reference, {best} has the largest estimated relative effect in this planning model. "
        "Interpret these values only after confirming extracted study data, risk of bias, transitivity, heterogeneity, and inconsistency."
    )
    return {"reference": reference, "direct": direct, "relative": relative, "league": league.reset_index(names="Treatment"), "narrative": narrative}


def dataframe_to_tiff(table: Any, title: str = "Table") -> bytes:
    """Render a dataframe as a TIFF image for download."""

    display_table = table.copy()
    try:
        import matplotlib.pyplot as plt

        fig_height = max(2.0, 0.35 * (len(display_table) + 2))
        fig_width = max(6.0, 1.4 * len(display_table.columns))
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=200)
        ax.axis("off")
        ax.set_title(title, fontsize=12, pad=10)
        mpl_table = ax.table(cellText=display_table.values, colLabels=display_table.columns, loc="center", cellLoc="center")
        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(8)
        mpl_table.scale(1, 1.25)
        buffer = BytesIO()
        fig.savefig(buffer, format="tiff", bbox_inches="tight")
        plt.close(fig)
        return buffer.getvalue()
    except ModuleNotFoundError:
        from PIL import Image, ImageDraw, ImageFont

        rows = [list(map(str, display_table.columns))] + display_table.astype(str).values.tolist()
        font = ImageFont.load_default()
        cell_w = 190
        cell_h = 28
        width = max(600, cell_w * max(1, len(display_table.columns)))
        height = cell_h * (len(rows) + 2)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)
        draw.text((8, 8), title, fill="black", font=font)
        y = cell_h
        for row_index, row in enumerate(rows):
            x = 0
            fill = "#EEF4FA" if row_index == 0 else "white"
            for cell in row:
                draw.rectangle([x, y, x + cell_w, y + cell_h], outline="black", fill=fill)
                draw.text((x + 5, y + 8), str(cell)[:28], fill="black", font=font)
                x += cell_w
            y += cell_h
        buffer = BytesIO()
        image.save(buffer, format="TIFF")
        return buffer.getvalue()


def plotly_figure_to_tiff(fig: Any, scale: int = 2) -> bytes:
    """Convert a Plotly figure to TIFF bytes via Kaleido PNG export and Pillow."""

    from PIL import Image

    png_bytes = fig.to_image(format="png", scale=scale)
    image = Image.open(BytesIO(png_bytes))
    output = BytesIO()
    image.save(output, format="TIFF")
    return output.getvalue()


def generate_protocol_notes(condition: str, interventions: Iterable[str], outcome: str) -> list[str]:
    """Return concise protocol reminders for a network meta-analysis workflow."""

    treatment_list = ", ".join(normalize_terms(interventions)) or "the eligible interventions"
    return [
        f"Define the population as patients with {condition or 'the target condition'} and document all eligibility criteria before screening.",
        f"Treatments to verify in full text: {treatment_list}. Combine doses/classes only if clinically justified a priori.",
        f"Primary outcome: {outcome or 'specify one primary outcome'}. Keep effect measures consistent across direct and indirect evidence.",
        "Use PubMed results as references for screening, but extract randomized study arms, sample sizes, and outcome data from full texts or trial reports.",
        "Assess transitivity, risk of bias, heterogeneity, inconsistency, and publication bias before interpreting rankings such as SUCRA or P-scores.",
    ]
