# PubMed Network Meta-analysis Agent

A Streamlit web application that helps prepare a network meta-analysis (NMA) by building a reproducible PICO-style PubMed query, retrieving PubMed-indexed papers through NCBI E-utilities, and drafting a treatment evidence network after manual screening.

> PubMed is an indexed biomedical literature database rather than a journal publisher. This project therefore uses papers indexed in PubMed as references and links each record to its official PubMed page.

## Features

- Builds transparent PubMed queries from condition, intervention, comparator, and outcome terms.
- Uses NCBI ESearch and EFetch endpoints to retrieve PubMed PMIDs, titles, journals, years, abstracts, publication types, and official PubMed URLs.
- Provides CSV export of retrieved references for screening.
- Flags treatment mentions in titles/abstracts as AI-assisted screening cues.
- Lets reviewers enter verified direct comparisons and visualizes the treatment network.
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
6. Export extracted study-arm outcome data to a validated NMA package such as `netmeta`, `gemtc`, BUGS/JAGS, or another pre-specified statistical workflow.

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
