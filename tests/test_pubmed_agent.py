import unittest

from pubmed_agent import (
    DISEASE_TEMPLATES,
    STUDY_DESIGNS,
    dataframe_to_tiff,
    build_pubmed_query,
    demo_extracted_comparisons,
    nma_dataset_from_edges,
    normalize_terms,
    pubmed_search_url,
    summarize_network,
)


class PubMedAgentTest(unittest.TestCase):
    def test_normalize_terms_accepts_common_separators(self):
        self.assertEqual(normalize_terms("adalimumab; etanercept\ninfliximab"), ["adalimumab", "etanercept", "infliximab"])

    def test_build_pubmed_query_contains_pico_and_filters(self):
        query = build_pubmed_query(
            condition="rheumatoid arthritis",
            interventions="methotrexate",
            comparators="placebo",
            outcome="ACR20 response",
        )

        self.assertIn('"rheumatoid arthritis"[Title/Abstract]', query)
        self.assertIn('"Arthritis, Rheumatoid"[MeSH Terms]', query)
        self.assertIn('"methotrexate"[Title/Abstract]', query)
        self.assertNotIn('"methotrexate"[MeSH Terms]', query)
        self.assertIn('"Placebos"[MeSH Terms]', query)
        self.assertIn('"ACR20 response"[Title/Abstract]', query)
        self.assertIn("randomized controlled trial[Publication Type]", query)
        self.assertIn("humans[MeSH Terms]", query)
        self.assertNotIn("network meta-analysis", query)

    def test_build_pubmed_query_prefers_title_abstract_for_unvalidated_condition_and_biologics(self):
        query = build_pubmed_query(
            condition="severe asthma",
            interventions="dupilumab, omalizumab, mepolizumab, benralizumab, tezepelumab, reslizumab",
            comparators="placebo/treatment",
            outcome="exacerbation",
            include_rct_filter=False,
            include_humans_english=False,
        )

        self.assertIn('"severe asthma"[Title/Abstract]', query)
        self.assertNotIn('"severe asthma"[MeSH Terms]', query)
        self.assertIn('"dupilumab"[Title/Abstract]', query)
        self.assertIn('"dupilumab"[Substance Name]', query)
        self.assertIn('"tezepelumab"[Substance Name]', query)
        self.assertIn('"placebo"[Title/Abstract]', query)
        self.assertNotIn("placebo/treatment", query)

    def test_asthma_covid_template_contains_required_terms_and_outcomes(self):
        template = DISEASE_TEMPLATES["Asthma and COVID-19"]
        query = build_pubmed_query(
            condition=template["condition"],
            interventions=template["interventions"],
            comparators=template["comparators"],
            outcome=template["outcome"],
            study_design="Observational Studies",
        )

        self.assertIn('"Asthma"[MeSH Terms]', query)
        self.assertIn('"COVID-19"[MeSH Terms]', query)
        self.assertIn('"ICS"[Title/Abstract]', query)
        self.assertIn('"ICS/LABA"[Title/Abstract]', query)
        self.assertIn('"oral corticosteroids"[Title/Abstract]', query)
        self.assertIn('"hospitalization"[Title/Abstract]', query)
        self.assertIn('"ICU admission"[Title/Abstract]', query)
        self.assertIn('"mechanical ventilation"[Title/Abstract]', query)
        self.assertIn('"mortality"[Title/Abstract]', query)
        self.assertIn("observational study[Publication Type]", query)
        self.assertNotIn("randomized controlled trial[Publication Type]", query)

    def test_study_design_both_includes_rct_and_observational_filters(self):
        query = build_pubmed_query("asthma", "dupilumab", "placebo", "exacerbation", study_design="Both")

        self.assertIn("randomized controlled trial[Publication Type]", query)
        self.assertIn("observational study[Publication Type]", query)
        self.assertIn("Both", STUDY_DESIGNS)

    def test_pubmed_search_url_points_to_pubmed(self):
        url = pubmed_search_url('"severe asthma"[Title/Abstract] AND "dupilumab"[Title/Abstract]')
        self.assertTrue(url.startswith("https://pubmed.ncbi.nlm.nih.gov/?term="))
        self.assertIn("%22severe+asthma%22%5BTitle%2FAbstract%5D", url)

    def test_summarize_network_returns_relative_effects_and_league(self):
        summary = summarize_network(demo_extracted_comparisons(), reference="placebo")

        self.assertEqual(summary["reference"], "placebo")
        self.assertFalse(summary["relative"].empty)
        self.assertIn("rank_score", summary["relative"].columns)
        self.assertIn("Treatment", summary["league"].columns)
        self.assertIn("planning model", summary["narrative"])

    def test_dataframe_to_tiff_returns_tiff_bytes(self):
        tiff = dataframe_to_tiff(demo_extracted_comparisons().head(1), title="Demo")

        self.assertTrue(tiff.startswith(b"II") or tiff.startswith(b"MM"))

    def test_nma_dataset_from_edges_has_export_columns(self):
        dataset = nma_dataset_from_edges(demo_extracted_comparisons().head(1))

        self.assertEqual(list(dataset.columns), ["PMID", "treatment", "comparator", "effect estimate", "standard error"])


if __name__ == "__main__":
    unittest.main()
