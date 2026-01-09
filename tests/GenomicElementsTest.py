
import unittest
import shutil
import os

import pandas as pd
import numpy as np

from .BedTableTest import TestBedTable3
from ..GenomicElements import GenomicElements
from ..BedTable import BedTable3

class TestGenomicElements(unittest.TestCase):
    def setUp(self):

        self.__wdir = "GenomicElement_test"

        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)

        self.__hg38_genome_path = "RGTools/large_files/hg38.fa"
        self.__bed3_region_file_path = os.path.join(self.__wdir, "test.bed3")

        TestBedTable3._gen_test_bed_file(self.__bed3_region_file_path)

        return super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.__wdir):
            shutil.rmtree(self.__wdir)
        return super().tearDown()

    def _init_GenomicElements(self):
        region_file_type = "bed3"
        region_file_path = self.__bed3_region_file_path
        genome_path = self.__hg38_genome_path
        return GenomicElements(region_file_path, region_file_type, genome_path)

    def test_one_hot_encoding(self):
        encoding = GenomicElements.one_hot_encoding("ACGT")
        self.assertTrue((encoding == np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])).all())
    
    def test_export_exogeneous_sequences(self):
        ge = self._init_GenomicElements()
        ge.export_exogeneous_sequences(os.path.join(self.__wdir, "test.fa"))

        with open(os.path.join(self.__wdir, "test.fa"), "r") as handle:
            lines = handle.readlines()
            self.assertEqual(len(lines), 8)
            self.assertEqual(lines[0], ">chr1:1-5\n")
            self.assertEqual(lines[1], "NNNN\n")
            self.assertEqual(lines[2], ">chr1:8-12\n")
            self.assertEqual(lines[3], "NNNN\n")

    def test_get_num_regions(self):
        ge = self._init_GenomicElements()
        self.assertEqual(ge.get_num_regions(), 4)

    def test_get_region_seq(self):
        ge = self._init_GenomicElements()
        bin1_seq = ge.get_region_seq("chr2", 127048023, 127107288)
        self.assertEqual(bin1_seq[:10], "ACGTTTTCTG")

        nonexist_seq = ge.get_region_seq("chr24", 127048023, 127107288)
        self.assertIsNone(nonexist_seq)
    
    def test_get_all_region_seqs(self):
        region_file_path = os.path.join(self.__wdir, "get_region_seq_test.bed3")
        region_file_type = "bed3"
        region_bt = BedTable3()
        region_bt.load_from_dataframe(pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [123400, 124400],
            "end": [123500, 124500], 
        }))
        region_bt.write(region_file_path)

        ge = GenomicElements(region_file_path, 
                             region_file_type, 
                             self.__hg38_genome_path, 
                             )

        region_seq_list = ge.get_all_region_seqs()
        self.assertEqual(len(region_seq_list), 2)
        self.assertEqual(region_seq_list[0][:12], "GTGTAATTACAA")

    def test_get_all_region_one_hot(self):
        region_file_path = os.path.join(self.__wdir, "one_hot_test.bed3")
        region_file_type = "bed3"
        region_bt = BedTable3()
        region_bt.load_from_dataframe(pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [123400, 124400],
            "end": [123500, 124500], 
        }))
        region_bt.write(region_file_path)

        ge = GenomicElements(region_file_path, 
                             region_file_type, 
                             self.__hg38_genome_path, 
                             )

        one_hot_arr = ge.get_all_region_one_hot()
        self.assertEqual(one_hot_arr.shape, (2, 100, 4))
    
    def test_load_region_anno_from_npy(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        npy_path = os.path.join(self.__wdir, "test.npy")
        anno_arr = np.array([1,2,3,4])
        np.save(npy_path, anno_arr)

        ge.load_region_anno_from_npy(anno_name, npy_path)
        self.assertEqual(ge.get_anno_dim(anno_name), 1)
        self.assertTrue((ge.get_anno_arr(anno_name).reshape(-1,) == anno_arr).all())

        anno_name = "test_anno_multi_dim"
        anno_arr = np.ones((4, 4))
        np.save(npy_path, anno_arr)
        ge.load_region_anno_from_npy(anno_name, npy_path)
        self.assertEqual(ge.get_anno_dim(anno_name), 4)
        self.assertTrue((ge.get_anno_arr(anno_name) == anno_arr).all())

        # Test loading from .npz file
        anno_name = "test_anno_npz"
        npz_path = os.path.join(self.__wdir, "test.npz")
        anno_arr = np.array([5, 6, 7, 8])
        np.savez(npz_path, anno_arr)
        
        ge.load_region_anno_from_npy(anno_name, npz_path)
        self.assertEqual(ge.get_anno_dim(anno_name), 1)
        self.assertTrue((ge.get_anno_arr(anno_name).reshape(-1,) == anno_arr).all())

        # Test error when .npz file contains multiple arrays
        multi_npz_path = os.path.join(self.__wdir, "test_multi.npz")
        arr1 = np.array([1, 2, 3, 4])
        arr2 = np.array([5, 6, 7, 8])
        np.savez(multi_npz_path, arr1=arr1, arr2=arr2)
        
        with self.assertRaises(ValueError) as context:
            ge.load_region_anno_from_npy("test_multi", multi_npz_path)
        self.assertIn("contains multiple arrays", str(context.exception))
        self.assertIn("Available keys", str(context.exception))
    
    def test_load_region_anno_from_arr(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])

        ge.load_region_anno_from_arr(anno_name, anno_arr)

        self.assertEqual(ge.get_anno_dim(anno_name), 1)
        self.assertTrue((ge.get_anno_arr(anno_name).reshape(-1,) == anno_arr).all())
        self.assertEqual(ge.get_anno_type(anno_name), "stat")

    def test_load_region_track_from_list(self):
        # Create a bed file with different lengths
        hetero_bed = os.path.join(self.__wdir, "hetero_load_list.bed")
        with open(hetero_bed, "w") as f:
            f.write("chr1\t123400\t123500\n") # 100bp
            f.write("chr1\t124400\t124450\n") # 50bp
        
        ge = GenomicElements(hetero_bed, "bed3", self.__hg38_genome_path)
        
        # Correct lengths: 100 and 50
        anno_list = [np.ones(100), np.ones(50)]
        ge.load_region_track_from_list("test_list", anno_list)
        
        output = ge.get_anno_arr("test_list")
        # Max length should be 100
        self.assertEqual(output.shape, (2, 100))
        # First region should be all ones
        self.assertTrue(np.all(output[0] == 1))
        # Second region should be all ones
        self.assertTrue(np.all(output[1, :50] == 1))
        self.assertTrue(np.all(output[1, 50:] == 0)) # padded area

        # Test get_region_anno_by_index
        self.assertTrue(np.all(ge.get_region_anno_by_index("test_list", 0) == np.ones(100)))
        self.assertTrue(np.all(ge.get_region_anno_by_index("test_list", 1) == np.ones(50)))

    def test_load_region_anno_from_arr_stat_shapes(self):
        ge = self._init_GenomicElements()
        # Accept (N,) and store as (N,1)
        anno = np.array([10, 20, 30, 40])
        ge.load_region_anno_from_arr("stat1", anno)
        self.assertEqual(ge.get_anno_arr("stat1").shape, (4, 1))
        self.assertEqual(ge.get_anno_type("stat1"), "stat")
        # Accept (N,1) unchanged
        anno2 = np.array([[1], [2], [3], [4]])
        ge.load_region_anno_from_arr("stat2", anno2)
        self.assertEqual(ge.get_anno_arr("stat2").shape, (4, 1))
        self.assertEqual(ge.get_anno_type("stat2"), "stat")

    def test_load_region_anno_from_arr_track_requires_max_len(self):
        ge = self._init_GenomicElements()
        # Track with correct width (max region length = 4)
        track = np.ones((4, 4))
        ge.load_region_anno_from_arr("track_ok", track)
        self.assertEqual(ge.get_anno_type("track_ok"), "track")
        self.assertEqual(ge.get_anno_dim("track_ok"), 4)
        # Track with wrong width should raise
        bad_track = np.ones((4, 3))
        with self.assertRaises(ValueError):
            ge.load_region_anno_from_arr("track_bad", bad_track)

    def test_get_region_anno_by_index(self):
        ge = self._init_GenomicElements()
        ge.load_region_anno_from_arr("stat", np.array([1, 2, 3, 4]))
        self.assertEqual(ge.get_region_anno_by_index("stat", 2), 3)
        ge.load_region_anno_from_arr("track", np.arange(16).reshape(4, 4))
        # Regions are length-homogeneous length=4 in test data
        self.assertTrue((ge.get_region_anno_by_index("track", 1) == np.array([4,5,6,7])).all())
    
    def test_save_anno_npy(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])
        ge.load_region_anno_from_arr(anno_name, anno_arr)

        npy_path = os.path.join(self.__wdir, "test.npy")
        ge.save_anno_npy(anno_name, npy_path)

        self.assertTrue(os.path.exists(npy_path))
        self.assertTrue((np.load(npy_path).reshape(-1,) == anno_arr).all())

    def test_apply_logical_filter(self):
        ge = self._init_GenomicElements()
        anno_name = "test_anno"
        anno_arr = np.array([1,2,3,4])
        ge.load_region_anno_from_arr(anno_name, anno_arr)

        # Track aligned to max len (4)
        track_arr = np.arange(16).reshape(4, 4)
        ge.load_region_anno_from_arr("track", track_arr)

        logical = np.array([True, False, True, False])

        new_region_file_path = os.path.join(self.__wdir, "new_region.bed3")
        new_ge = ge.apply_logical_filter(logical, new_region_file_path)

        self.assertEqual(new_ge.get_num_regions(), 2)
        self.assertEqual(new_ge.get_anno_dim(anno_name), 1)
        self.assertTrue((new_ge.get_anno_arr(anno_name).reshape(-1,) == np.array([1,3])).all())
        self.assertEqual(new_ge.get_anno_type(anno_name), "stat")

        # Track should be trimmed to new max len after filter (regions kept are same length=4 here)
        self.assertEqual(new_ge.get_anno_type("track"), "track")
        self.assertEqual(new_ge.get_anno_arr("track").shape, (2, 4))
        np.testing.assert_array_equal(new_ge.get_anno_arr("track")[0], np.array([0,1,2,3]))
        np.testing.assert_array_equal(new_ge.get_anno_arr("track")[1], np.array([8,9,10,11]))

