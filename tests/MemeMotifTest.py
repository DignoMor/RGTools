
import unittest
import requests
import shutil
import os

from ..MemeMotif import MemeMotif

class TestMemeMotif(unittest.TestCase):
    def setUp(self):
        self._test_dir = "MemeMotif_test"

        if not os.path.exists(self._test_dir):
            os.makedirs(self._test_dir)

        self._meme_file_path = os.path.join(self._test_dir, "test_motifs.meme")
        # Download meme example from meme website
        import requests

        url = "https://meme-suite.org/meme/doc/examples/sample-dna-motif.meme"
        response = requests.get(url)

        # Save the file locally
        with open(self._meme_file_path, "wb") as f:
            f.write(response.content)


    def tearDown(self):
        if os.path.exists(self._test_dir):
            shutil.rmtree(self._test_dir)

        return super().tearDown()

    def test_get_meme_version(self):
        meme = MemeMotif(self._meme_file_path)
        version = meme.get_meme_version()
        self.assertEqual(version, "4")
    
    def test_get_alphabet(self):
        meme = MemeMotif(self._meme_file_path)
        alphabet = meme.get_alphabet()
        self.assertEqual(alphabet, "ACGT")
    
    def test_get_strands(self):
        meme = MemeMotif(self._meme_file_path)
        strands = meme.get_strands()
        self.assertEqual(len(strands), 2)
        self.assertEqual(strands[0], "+")
        self.assertEqual(strands[1], "-")

    def test_get_bg_freq(self):
        meme = MemeMotif(self._meme_file_path)
        bg_freq = meme.get_bg_freq()
        self.assertEqual(len(bg_freq), 4)
        self.assertAlmostEqual(bg_freq[0], 0.303)
        self.assertAlmostEqual(bg_freq[1], 0.183)
        self.assertAlmostEqual(bg_freq[2], 0.209)
        self.assertAlmostEqual(bg_freq[3], 0.306)

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
