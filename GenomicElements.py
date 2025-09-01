
import os 

import numpy as np

from Bio import SeqIO

from .BedTable import BedTable3, BedTable6, BedTable6Plus, BedTable3Plus
from .GeneralElements import GeneralElements

class GenomicElements(GeneralElements):
    '''
    Class for Genomics elements.
    '''

    def __init__(self, region_path, region_file_type, fasta_path):
        '''
        Constructor for GenomicElements class.

        Keyword arguments:
        - region_path: Path to the region file.
        - region_file_type: Type of the region file (see GenomicElements.get_region_file_suffix2class_dict).
        - fasta_path: Path to the genome file.
        '''
        super().__init__()
        self.region_path = region_path
        self._region_file_type = region_file_type
        if not self._region_file_type in self.get_region_file_suffix2class_dict().keys():
            raise ValueError(f"Invalid region file type: {self._region_file_type}")

        self._fasta_path = fasta_path

    @property
    def fasta_path(self):
        return self._fasta_path

    @property
    def region_file_type(self):
        return self._region_file_type

    @staticmethod
    def get_region_file_suffix2class_dict():
        '''
        Return the dictionary that maps region file suffix to the corresponding class.
        '''
        return {
            "bed3": BedTable3,
            "bed6": BedTable6,
            "bed6gene": GenomicElements.BedTable6Gene,
            "bed3gene": GenomicElements.BedTable3Gene,
        }

    @staticmethod
    def BedTable6Gene(enable_sort=True):
        '''
        Helper function to return a BedTable6Plus object
        that can load bed6gene annotations.
        '''
        bt = BedTable6Plus(extra_column_names=["gene_symbol"], 
                           extra_column_dtype=[str], 
                           enable_sort=enable_sort,
                           )
        return bt

    @staticmethod
    def BedTable3Gene(enable_sort=True):
        '''
        Helper function to return a BedTable3Plus object
        that can load bed3gene annotations.
        '''
        bt = BedTable3Plus(extra_column_names=["gene_symbol"], 
                           extra_column_dtype=[str], 
                           enable_sort=enable_sort,
                           )
        return bt
    
    @staticmethod
    def set_parser_genome(parser):
        parser.add_argument("--fasta_path", 
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

    def export_exogeneous_sequences(self, fasta_path):
        '''
        Export regions as exogeneous sequences.

        Keyword arguments:
        - fasta_path: Path to save the exogeneous sequences.

        Returns:
        - None
        '''
        if os.path.exists(fasta_path):
            raise ValueError(f"File {fasta_path} already exists.")
        
        for region in self.get_region_bed_table().iter_regions():
            seq = self.get_region_seq(region["chrom"], region["start"], region["end"])
            with open(fasta_path, "a") as handle:
                handle.write(f">{region['chrom']}:{region['start']}-{region['end']}\n{seq}\n")

    def get_region_bed_table(self):
        '''
        Return a bed table object for the region file.
        '''
        bt = self.get_region_file_suffix2class_dict()[self.region_file_type](enable_sort=False)
        bt.load_from_file(self.region_path)

        return bt
    
    def get_region_seq(self, chrom: str, start: int, end: int) -> str:
        '''
        Return the sequene of a given region. 
        The coordinates are of bed convention.
        
        Keyword arguments:
        - chrom: Chromosome name.
        - start: Start coordinate.
        - end: End coordinate.
        
        Returns:
        - seq: Sequence of the region.
        '''
        with open(self.fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == chrom:
                    return str(record.seq[start:end])  # Convert to 0-based index
        return None

    def apply_logical_filter(self, logical, new_region_path):
        '''
        Apply logical filter to the regions.

        Keyword arguments:
        - logical: np.Array, Logical array to filter the regions.
        - new_region_path: Path to save the new region_file for filtered regions.

        Returns:
        - a new GenomicElements object with the filtered regions.
        '''
        result_ge = self.__class__(region_path=new_region_path,
                                   region_file_type=self.region_file_type,
                                   fasta_path=self.fasta_path,
                                   )
        
        new_bt = self.get_region_bed_table().apply_logical_filter(logical)
        new_bt.write(new_region_path)

        for anno_name, anno_arr in self._anno_arr_dict.items():
            new_anno_arr = anno_arr[logical]
            result_ge.load_region_anno_from_arr(anno_name, new_anno_arr)

        return result_ge
    