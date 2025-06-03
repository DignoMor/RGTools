
import unittest
import shutil
import os

from RGTools.ExogeneousSequences import ExogeneousSequences
from RGTools.GenomicElements import GenomicElements

class TestExogeneousSequences(unittest.TestCase):
    def setUp(self):
        self.__wdir = "ExogeneousSequences_test"

        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)

        self.__fasta_path = os.path.join(self.__wdir, "test.fa")

        with open(self.__fasta_path, "w") as handle:
            handle.write(">chr2:127048023-127048033\nACGTTTTCTG\n")
            handle.write(">chr1:123500-123512\nGTGTAATTACAA\n")

        return super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.__wdir):
            shutil.rmtree(self.__wdir)

    def test_get_region_bed_table(self):
        es = ExogeneousSequences(self.__fasta_path)
        bt = es.get_region_bed_table()
        bt.write(os.path.join(self.__wdir, "test.bed3"))

        ge = GenomicElements(region_path=os.path.join(self.__wdir, "test.bed3"),
                             region_file_type="bed3", 
                             fasta_path=self.__fasta_path, 
                             )
        seqs = []
        for region in bt.iter_regions():
            seqs.append(ge.get_region_seq(region["chrom"], region["start"], region["end"]))
        
        self.assertEqual(seqs[0], "ACGTTTTCTG")
