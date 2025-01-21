
import unittest
import shutil
import os

import numpy as np

from ..BwTrack import BwTrack

class TestBwTrack(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    
    def tearDown(self):
        return super().tearDown()
    
    def test_quantify_signal(self):
        input_signal = np.array([1,2,3,4,5])

        self.assertEqual(BwTrack.quantify_signal(input_signal, "raw_count"), 
                         input_signal.sum(), 
                         )
        
        self.assertEqual(BwTrack.quantify_signal(input_signal, "RPK"),
                         input_signal.sum() / len(input_signal) * 1e3,
                         )

