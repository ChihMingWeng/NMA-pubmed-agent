import unittest

from pubmed_agent import (
    dataframe_to_tiff,
    build_pubmed_query,
    demo_extracted_comparisons,
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
            interventions="adalimumab, etanercept",
            comparators="placebo",
            outcome="ACR20 response",
        )

        self.assertIn('"rheumatoid arthritis"[MeSH Terms]', query)
        self.assertIn('"adalimumab"[Title/Abstract]', query)
        self.assertIn('"placebo"[MeSH Terms]', query)
        self.assertIn('"ACR20 response"[Title/Abstract]', query)
        self.assertIn("randomized controlled trial[Publication Type]", query)
        self.assertIn("humans[MeSH Terms]", query)

    def test_pubmed_search_url_points_to_pubmed(self):
        self.assertTrue(pubmed_search_url("asthma").startswith("https://pubmed.ncbi.nlm.nih.gov/?term="))

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


if __name__ == "__main__":
    unittest.main()
