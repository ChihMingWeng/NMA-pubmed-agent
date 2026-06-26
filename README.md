# PubMed Network Meta-analysis Agent

A Streamlit web application that helps prepare a network meta-analysis (NMA) by building a reproducible PICO-style PubMed query, retrieving PubMed-indexed papers through NCBI E-utilities, and drafting a treatment evidence network after manual screening.

> PubMed is an indexed biomedical literature database rather than a journal publisher. This project therefore uses papers indexed in PubMed as references and links each record to its official PubMed page.

## Features

- Builds transparent PubMed queries from condition, intervention, comparator, and outcome terms.
- Provides disease templates for rheumatoid arthritis, severe asthma biologics, and asthma with COVID-19.
- Supports randomized controlled trial, observational study, or combined study-design search filters.
- Uses NCBI ESearch and EFetch endpoints to retrieve PubMed PMIDs, titles, journals, years, abstracts, publication types, and official PubMed URLs.
- Provides CSV export of retrieved references for screening.
- Flags treatment mentions in titles/abstracts as AI-assisted screening cues.
- Lets reviewers enter verified direct comparisons and visualizes the treatment network.
- Ships with a runnable severe-asthma biologics PubMed topic and synthetic extracted contrasts so the web workflow can be explored immediately.
- Produces an exploratory relative-effect summary, forest plot, league table, and narrative interpretation from entered log-effect estimates and standard errors.
- Exports an NMA-ready dataset table with PMID, treatment, comparator, effect estimate, and standard error columns.
- Exports the network plot, forest plot, and league table as TIFF (`.tif`) images for manuscript or slide preparation.
- Generates protocol reminders for transitivity, risk of bias, heterogeneity, inconsistency, and interpretation of treatment rankings.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Recommended workflow

1. Enter the clinical condition, interventions, comparators, and outcome.
2. Review and edit the generated PubMed query if needed.
3. Search PubMed and download the references CSV.
4. Screen titles/abstracts and full texts using at least two reviewers where appropriate.
5. Enter only verified direct treatment comparisons into the network editor.
6. Review the exploratory relative-effect table, forest plot, league table, and generated text.
7. Download TIFF versions of the figures/tables if needed.
8. Export extracted study-arm outcome data to a validated NMA package such as `netmeta`, `gemtc`, BUGS/JAGS, or another pre-specified statistical workflow.

## Statistical note

The built-in summary pools duplicate direct contrasts with inverse-variance weighting and combines simple network paths to give a transparent planning approximation. It is designed for teaching, protocol development, and data-checking only. Final clinical conclusions should use a validated frequentist or Bayesian NMA implementation with prespecified effect measures, heterogeneity, inconsistency, and risk-of-bias handling.

## NCBI usage

NCBI recommends including an email address with E-utilities requests. The app has an optional email field that is passed to NCBI when provided.

## Tests and checks

```bash
python -m compileall app.py pubmed_agent.py
python - <<'PY'
from pubmed_agent import build_pubmed_query, pubmed_search_url
query = build_pubmed_query('asthma', 'budesonide, fluticasone', 'placebo', 'exacerbation')
assert 'asthma' in query
assert 'budesonide' in query
assert pubmed_search_url(query).startswith('https://pubmed.ncbi.nlm.nih.gov/?term=')
print(query)
PY
```
