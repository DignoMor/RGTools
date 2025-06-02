
import unittest
import requests
import shutil
import os

import numpy as np

from ..MemeMotif import MemeMotif

class TestMemeMotif(unittest.TestCase):
    def setUp(self):
        self._test_dir = "MemeMotif_test"

        if not os.path.exists(self._test_dir):
            os.makedirs(self._test_dir)

        self._meme_file_path = os.path.join(self._test_dir, "test_motifs.meme")
        # Download meme example from meme website
        url = "https://meme-suite.org/meme/doc/examples/sample-dna-motif.meme"
        response = requests.get(url)

        # Save the file locally
        with open(self._meme_file_path, "wb") as f:
            f.write(response.content)

    def tearDown(self):
        if os.path.exists(self._test_dir):
            shutil.rmtree(self._test_dir)

        return super().tearDown()
    
    def test_write_meme_file(self):
        meme = MemeMotif(self._meme_file_path)
        meme.write_meme_file(os.path.join(self._test_dir, "test_motifs_copy.meme"))

        output_meme = MemeMotif(os.path.join(self._test_dir, "test_motifs_copy.meme"))
        self.assertEqual(meme.get_meme_version(), output_meme.get_meme_version())
        self.assertEqual(meme.get_alphabet(), output_meme.get_alphabet())
        self.assertEqual(meme.get_strands(), output_meme.get_strands())
        self.assertEqual(meme.get_bg_freq(), output_meme.get_bg_freq())
        self.assertEqual(meme.get_motif_list(), output_meme.get_motif_list())
        self.assertTrue((meme.get_motif_pwm("crp") == output_meme.get_motif_pwm("crp")).all())
        self.assertEqual(meme.get_motif_alphabet_length("crp"), output_meme.get_motif_alphabet_length("crp"))
        self.assertEqual(meme.get_motif_length("crp"), output_meme.get_motif_length("crp"))
        self.assertEqual(meme.get_motif_num_source_sites("crp"), output_meme.get_motif_num_source_sites("crp"))
        self.assertEqual(meme.get_motif_source_eval("crp"), output_meme.get_motif_source_eval("crp"))

    def test_get_meme_version(self):
        meme = MemeMotif(self._meme_file_path)
        version = meme.get_meme_version()
        self.assertEqual(version, "4")
    
    def test_set_meme_version(self):
        meme = MemeMotif(self._meme_file_path)
        meme.set_meme_version("5")
        self.assertEqual(meme.get_meme_version(), "5")

    def test_get_alphabet(self):
        meme = MemeMotif(self._meme_file_path)
        alphabet = meme.get_alphabet()
        self.assertEqual(alphabet, "ACGT")
    
    def test_set_alphabet(self):
        meme = MemeMotif(self._meme_file_path)
        meme.set_alphabet("ACGTN")
        self.assertEqual(meme.get_alphabet(), "ACGTN")

    def test_get_strands(self):
        meme = MemeMotif(self._meme_file_path)
        strands = meme.get_strands()
        self.assertEqual(len(strands), 2)
        self.assertEqual(strands[0], "+")
        self.assertEqual(strands[1], "-")

    def test_set_strands(self):
        meme = MemeMotif(self._meme_file_path)
        meme.set_strands(["+", "-", "."])
        self.assertEqual(meme.get_strands(), ["+", "-", "."])

    def test_get_bg_freq(self):
        meme = MemeMotif(self._meme_file_path)
        bg_freq = meme.get_bg_freq()
        self.assertEqual(len(bg_freq), 4)
        self.assertAlmostEqual(bg_freq[0], 0.303)
        self.assertAlmostEqual(bg_freq[1], 0.183)
        self.assertAlmostEqual(bg_freq[2], 0.209)
        self.assertAlmostEqual(bg_freq[3], 0.306)
    
    def test_set_bg_freq(self):
        meme = MemeMotif(self._meme_file_path)
        meme.set_bg_freq({"A": 0.25, "C": 0.25, "G": 0.15, "T": 0.35})
        self.assertEqual(meme.get_bg_freq(), {"A": 0.25, "C": 0.25, "G": 0.15, "T": 0.35})

    def test_get_motif_list(self):
        meme = MemeMotif(self._meme_file_path)
        motifs = meme.get_motif_list()
        self.assertEqual(len(motifs), 2)
        self.assertEqual(motifs[0], "crp")
        self.assertEqual(motifs[1], "lexA")
    
    def test_get_motif_pwm(self):
        meme = MemeMotif(self._meme_file_path)
        motif_pwm = meme.get_motif_pwm("crp")
        self.assertEqual(motif_pwm.shape, (19, 4))
        self.assertEqual(motif_pwm[0, 0], 0)
        self.assertEqual(motif_pwm[3, 2], 0.764706)
    
    def test_get_motif_alphabet_length(self):
        meme = MemeMotif(self._meme_file_path)

        alphabet_length = meme.get_motif_alphabet_length("lexA")
        self.assertEqual(alphabet_length, 4)
    
    def test_get_motif_length(self):
        meme = MemeMotif(self._meme_file_path)
        motif_length = meme.get_motif_length("lexA")
        self.assertEqual(motif_length, 18)
    
    def test_get_motif_num_source_sites(self):
        meme = MemeMotif(self._meme_file_path)
        motif_num = meme.get_motif_num_source_sites("lexA")
        self.assertEqual(motif_num, 14)
    
    def test_get_motif_source_eval(self):
        meme = MemeMotif(self._meme_file_path)
        source_eval = meme.get_motif_source_eval("crp")
        self.assertEqual(source_eval, 4.1e-9)
    
    def test_add_motif(self):
        meme = MemeMotif(self._meme_file_path)
        dpe_pwm = np.array([[0.46808, 0.06382, 0.383, 0.0851],
                            [0.0, 0.21, 0.79, 0.0],
                            [0.5745, 0.0425, 0.021, 0.362],
                            [0.043, 0.489, 0.021, 0.447],
                            [0.064, 0.255, 0.66, 0.021],
                            [0.064, 0.1702, 0.1915, 0.5743]])

        meme.add_motif("dpe", {"alphabet_length": 4, 
                               "motif_length": 6, 
                               "num_source_sites": 1000, 
                               "source_eval": 1, 
                               "pwm": dpe_pwm})

        self.assertEqual(meme.get_motif_list(), ["crp", "lexA", "dpe"])
        self.assertEqual(meme.get_motif_alphabet_length("dpe"), 4)
        self.assertEqual(meme.get_motif_length("dpe"), 6)
        self.assertEqual(meme.get_motif_num_source_sites("dpe"), 1000)
        self.assertEqual(meme.get_motif_source_eval("dpe"), 1)
        self.assertTrue((meme.get_motif_pwm("dpe") == dpe_pwm).all())

    def test_calculate_pwm_score(self):
        meme = MemeMotif(self._meme_file_path)
        motif_pwm = meme.get_motif_pwm("crp")
        score = MemeMotif.calculate_pwm_score("CGCGATCGATCGTTAAGTT", 
                                              motif_pwm, 
                                              bg_freq=meme.get_bg_freq(), 
                                              )
        self.assertAlmostEqual(score, -0.45161817)
        score = MemeMotif.calculate_pwm_score("CGCGATCGATCGTTAAGTT", 
                                              motif_pwm, 
                                              bg_freq=meme.get_bg_freq(), 
                                              reverse_complement=True, 
                                              )
        self.assertAlmostEqual(score, -np.inf)

    def test_search_one_motif(self):
        meme = MemeMotif(self._meme_file_path)
        score_arr = MemeMotif.search_one_motif("CGCGATCGATCGTTAAGTTG", 
                                               meme.get_alphabet(), 
                                               meme.get_motif_pwm("crp"), 
                                               bg_freq=meme.get_bg_freq(), 
                                               )
        self.assertEqual(len(score_arr), 20)
        self.assertAlmostEqual(score_arr[1], -np.inf)
        self.assertAlmostEqual(score_arr[0], -0.45161817)
