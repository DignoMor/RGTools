
import numpy as np

from .BedTable import BedTable3, BedTable6, BedTable6Plus

class GenomicElements:
    '''
    Class for Genomics elements.
    '''

    def __init__(self, region_path, region_file_type, genome_path):
        '''
        Constructor for GenomicElements class.

        Keyword arguments:
        - region_path: Path to the region file.
        - region_file_type: Type of the region file (see GenomicElements.get_region_file_suffix2class_dict).
        - genome_path: Path to the genome file.
        '''
        self.region_path = region_path
        self.region_file_type = region_file_type
        if not self.region_file_type in self.get_region_file_suffix2class_dict().keys():
            raise ValueError(f"Invalid region file type: {self.region_file_type}")

        self._anno_arr_dict = {}
        self._anno_length_dict = {}
        self.genome_path = genome_path

    @staticmethod
    def get_region_file_suffix2class_dict():
        '''
        Return the dictionary that maps region file suffix to the corresponding class.
        '''
        return {
            "bed3": BedTable3,
            "bed6": BedTable6,
            "bed6gene": GenomicElements.BedTable6Gene,
        }

    @staticmethod
    def BedTable6Gene():
        '''
        Helper function to return a BedTable6Plus object
        that can load bed6gene annotations.
        '''
        bt = BedTable6Plus(extra_column_names=["gene_symbol"], 
                        extra_column_dtype=[str], 
                        )
        return bt
    
    @staticmethod
    def set_parser_genome(parser):
        parser.add_argument("--genome_path", 
                            help="Path to the genome file.",
                            required=True,
                            type=str, 
                            )

    @staticmethod
    def set_parser_genomic_element_region(parser):
        parser.add_argument("--region_file_path", 
                            help="Path to the region file.",
                            required=True,
                            type=str, 
                            )
        
        parser.add_argument("--region_file_type",
                            help="Type of the region file. "
                                 "Valid types: {}".format(
                                     list(GenomicElements.get_region_file_suffix2class_dict().keys())
                                     ),
                            required=True,
                            default="bed3", 
                            type=str, 
                            choices=GenomicElements.get_region_file_suffix2class_dict().keys(),
                            )

    def get_region_bed_table(self):
        '''
        Return a bed table object for the region file.
        '''
        bt = self.get_region_file_suffix2class_dict()[self.region_file_type](enable_sort=False)
        bt.load_from_file(self.region_path)

        return bt
    
    def get_num_regions(self):
        '''
        Return the number of regions.
        '''
        return self.get_region_bed_table().__len__()
    
    def load_region_anno_from_npy(self, anno_name, npy_path):
        '''
        Load annotation for each element in the region file. 
        
        Keyword arguments:
        - anno_name: Name of the annotation.
        - npy_path: Path to the numpy file.

        Returns:
        None
        '''
        anno_arr = np.load(npy_path)
        self.load_region_anno_from_arr(anno_name, anno_arr)

    def load_region_anno_from_arr(self, anno_name, anno_arr):
        '''
        Load annotation for each element in the region file. 
        
        Keyword arguments:
        - anno_name: Name of the annotation.
        - anno_arr: np.Array, Annotation array.

        Returns:
        None
        '''
        anno_shape = anno_arr.shape

        if anno_shape[0] != self.get_num_regions():
            raise ValueError(f"Annotation shape {anno_shape} does not match the number of regions: {self.get_num_regions()}")
        
        if len(anno_shape) == 1:
            self._anno_arr_dict[anno_name] = anno_arr
            self._anno_length_dict[anno_name] = 1
        
        else:
            self._anno_arr_dict[anno_name] = anno_arr
            self._anno_length_dict[anno_name] = anno_shape[1]

    def get_anno_dim(self, anno_name):
        '''
        Return the dimension of the annotation.

        Keyword arguments:
        - anno_name: Name of the annotation.

        Return: 
        - dim: Dimension of the annotation.
        '''
        return self._anno_length_dict[anno_name]

    def get_anno_arr(self, anno_name):
        '''
        Return the annotation array.

        Keyword arguments:
        - anno_name: Name of the annotation.

        Return:
        - anno_arr: np.Array, Annotation array.
        '''
        return self._anno_arr_dict[anno_name]
    
    def apply_logical_filter(self, logical, new_region_path):
        '''
        Apply logical filter to the regions.

        Keyword arguments:
        - logical: np.Array, Logical array to filter the regions.
        - region_path: Path to save the new region_file for filtered regions.

        Returns:
        - a new GenomicElements object with the filtered regions.
        '''
        result_ge = self.__class__(region_path=new_region_path,
                                   region_file_type=self.region_file_type,
                                   genome_path=self.genome_path,
                                   )
        
        new_bt = self.get_region_bed_table().apply_logical_filter(logical)
        new_bt.write(new_region_path)

        for anno_name, anno_arr in self._anno_arr_dict.items():
            new_anno_arr = anno_arr[logical]
            result_ge.load_region_anno_from_arr(anno_name, new_anno_arr)

        return result_ge
    