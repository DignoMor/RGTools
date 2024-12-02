
import unittest

from ..SNP_utils import EnsemblRestSearch

class TestEnsemblRestSearch(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    
    def tearDown(self):
        return super().tearDown()

    def test_get_rsid_info(self):

        rsid = "rs56116432"

        search_engine = EnsemblRestSearch(genome_version="hg38")
        info = search_engine.get_rsid_info(rsid)

        self.assertEqual(info['name'], rsid)
        self.assertEqual(info["mappings"][0]["location"], "9:133256042-133256042")
        self.assertEqual(info["mappings"][0]["assembly_name"], "GRCh38")
        self.assertEqual(info["mappings"][0]["seq_region_name"], "9")
        self.assertEqual(info["mappings"][0]["start"], 133256042)

    def test_get_rsid_from_location(self):

        chrom = "1"
        pos = 161155392

        search_engine = EnsemblRestSearch(genome_version="hg19")
        rsids = search_engine.get_rsid_from_location(chrom, pos)

        self.assertEqual(len(rsids), 1)
        self.assertIn("rs4575098", rsids)

        chrom = "chr1"
        pos = 630947
        search_engine = EnsemblRestSearch(genome_version="hg38")
        rsids = search_engine.get_rsid_from_location(chrom, pos)

        self.assertEqual(len(rsids), 2)
        self.assertIn("rs1390538076", rsids)
        self.assertIn("rs2100353415", rsids)

        chrom = "chr1"
        pos = 17957656
        search_engine = EnsemblRestSearch(genome_version="hg38")
        rsids = search_engine.get_rsid_from_location(chrom, pos)

        self.assertEqual(len(rsids), 1)

    def test_is_SNP(self):
        rsid = "rs1367061710"
        search_engine = EnsemblRestSearch(genome_version="hg38")
        is_snp, simple_snp_info = search_engine._is_SNP(rsid)
        self.assertFalse(is_snp)

        rsid = "rs4575098"
        search_engine = EnsemblRestSearch(genome_version="hg38")
        is_snp, simple_snp_info = search_engine._is_SNP(rsid)
        self.assertTrue(is_snp)

    def test_get_rsid_snp_simple_info(self):

        rsid = "rs56116432"

        search_engine = EnsemblRestSearch(genome_version="hg38")
        info = search_engine.get_rsid_snp_simple_info(rsid)

        self.assertEqual(info["chrom"], "chr9")
        self.assertEqual(info["start"], 133256042)
        self.assertEqual(info["end"], 133256043)
        self.assertEqual(info["bases"], "C/A/T")
    
    def test_prioritize_rsids(self):
        rsids = ["rs1390538076", "rs2100353415"]

        search_engine = EnsemblRestSearch(genome_version="hg38")
        prioritized_rsid, prioritized_snp_info = search_engine.prioritize_rsids(rsids)

        self.assertEqual(prioritized_rsid, "rs1390538076")
        self.assertEqual(prioritized_snp_info["chrom"], "chr1")
        self.assertEqual(prioritized_snp_info["start"], 630947)
        self.assertEqual(prioritized_snp_info["end"], 630948)

        rsids = ["rs1367061710", "rs4575098"]
        prioritized_rsid, prioritized_snp_info = search_engine.prioritize_rsids(rsids)

        self.assertEqual(prioritized_rsid, "rs4575098")
        self.assertEqual(prioritized_snp_info["chrom"], "chr1")
        self.assertEqual(prioritized_snp_info["start"], 161185602)
        self.assertEqual(prioritized_snp_info["end"], 161185603)
    
