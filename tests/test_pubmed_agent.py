import unittest

from pubmed_agent import build_pubmed_query, normalize_terms, pubmed_search_url


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


if __name__ == "__main__":
    unittest.main()
