
import unittest
import os

import numpy as np

from .BedTableTest import TestBedTable3
from ..GenomicElements import GenomicElements

class TestGenomicElements(unittest.TestCase):
    def setUp(self):

        self.__wdir = "GenomicElement_test"

        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)

        self.__hg38_genome_path = "large_files/hg38.fa"
        self.__bed3_region_path = os.path.join(self.__wdir, "test.bed3")

        TestBedTable3._gen_test_bed_file(self.__bed3_region_path)

        return super().setUp()
    
    def tearDown(self):
        return super().tearDown()

    def _init_GenomicElements(self):
        region_path = self.__bed3_region_path
        region_file_type = "bed3"
        genome_path = self.__hg38_genome_path
        return GenomicElements(region_path, region_file_type, genome_path)

    def test_get_num_regions(self):
        ge = self._init_GenomicElements()
        self.assertEqual(ge.get_num_regions(), 4)
    
    def test_load_region_anno_from_npy(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        npy_path = os.path.join(self.__wdir, "test.npy")
        anno_arr = np.array([1,2,3,4])
        np.save(npy_path, anno_arr)

        ge.load_region_anno_from_npy(anno_name, npy_path)
        self.assertEqual(ge.get_anno_dim(anno_name), 1)
        self.assertTrue((ge.get_anno_arr(anno_name) == anno_arr).all())

        anno_name = "test_anno_multi_dim"
        anno_arr = np.ones((4, 4))
        np.save(npy_path, anno_arr)
        ge.load_region_anno_from_npy(anno_name, npy_path)
        self.assertEqual(ge.get_anno_dim(anno_name), 4)
        self.assertTrue((ge.get_anno_arr(anno_name) == anno_arr).all())
    
    def test_load_region_anno_from_arr(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])

        ge.load_region_anno_from_arr(anno_name, anno_arr)

        self.assertEqual(ge.get_anno_dim(anno_name), 1)
        self.assertTrue((ge.get_anno_arr(anno_name) == anno_arr).all())
    
    def test_save_anno_npy(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])
        ge.load_region_anno_from_arr(anno_name, anno_arr)

        npy_path = os.path.join(self.__wdir, "test.npy")
        ge.save_anno_npy(anno_name, npy_path)

        self.assertTrue(os.path.exists(npy_path))
        self.assertTrue((np.load(npy_path) == anno_arr).all())

    def test_apply_logical_filter(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])
        ge.load_region_anno_from_arr(anno_name, anno_arr)

        logical = np.array([True, False, True, False])

        new_region_path = os.path.join(self.__wdir, "new_region.bed3")
        new_ge = ge.apply_logical_filter(logical, new_region_path)

        self.assertEqual(new_ge.get_num_regions(), 2)
        self.assertEqual(new_ge.get_anno_dim(anno_name), 1)
        self.assertTrue((new_ge.get_anno_arr(anno_name) == np.array([1,3])).all())

