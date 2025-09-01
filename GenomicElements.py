
import os 
import abc
import numpy as np

from Bio import SeqIO

from .BedTable import BedTable3, BedTable6, BedTable6Plus, BedTable3Plus

class GeneralElements(abc.ABC):
    '''
    Abstract base class for GenomicElements and ExogeneousSequences.
    
    This class provides common functionality for handling genomic elements
    and exogeneous sequences, including:
    - Annotation management
    - Sequence operations
    - Bed table operations
    - Common utility methods
    '''

    def __init__(self):
        '''
        Initialize common attributes for all element types.
        '''
        self._anno_arr_dict = {}
        self._anno_length_dict = {}

    @property
    @abc.abstractmethod
    def fasta_path(self):
        '''
        Abstract property for the path to the fasta file.
        Must be implemented by subclasses.
        '''
        pass

    @property
    @abc.abstractmethod
    def region_file_type(self):
        '''
        Abstract property for the region file type.
        Must be implemented by subclasses.
        '''
        pass

    @abc.abstractmethod
    def get_region_bed_table(self):
        '''
        Abstract method to return a bed table object for the regions.
        Must be implemented by subclasses.
        '''
        pass

    def get_all_region_seqs(self):
        '''
        Get the sequences for all regions.
        This function reads the genome file to memory.
        
        Returns:
        - A list of strings of length region_length
        '''
        # get chrom name to seq dict
        chrom_seq_dict = {}
        with open(self.fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                chrom_seq_dict[record.id] = str(record.seq)

        out_seqs = []

        for i, region in enumerate(self.get_region_bed_table().iter_regions()):
            if not region["chrom"] in chrom_seq_dict:
                raise ValueError(f"Chromosome {region['chrom']} not found in the genome file.")
            seq = chrom_seq_dict[region["chrom"]][region["start"]:region["end"]]
            out_seqs.append(seq)
        
        return out_seqs

    def get_all_region_one_hot(self):
        '''
        Get the one hot encoding for all regions.
        This function reads the genome file to memory 
        so that it is memory intensive.
        '''
        out_seqs = self.get_all_region_seqs()

        out_seq_lens = np.array([len(s) for s in out_seqs])
        if not (out_seq_lens == out_seq_lens[0]).all():
            raise ValueError(f"Region sequences have different lengths: {out_seq_lens}")
        
        out_arr = np.zeros((len(out_seqs), out_seq_lens[0], 4), dtype="int8")
        for i, seq in enumerate(out_seqs):
            out_arr[i] = self.one_hot_encoding(seq)
        return out_arr

    @abc.abstractmethod
    def apply_logical_filter(self, logical, new_path):
        '''
        Abstract method to apply logical filter to the regions.
        Must be implemented by subclasses.
        '''
        pass

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
    
    def save_anno_npy(self, anno_name, npy_path):
        '''
        Save annotation to a numpy file.

        Keyword arguments:
        - anno_name: Name of the annotation.
        - npy_path: Path to save the annotation.

        Returns:
        None
        '''
        np.save(npy_path, self.get_anno_arr(anno_name))

    @staticmethod
    def one_hot_encoding(seq: str):
        '''
        Return the one hot encoding of a sequence.
        Ambiguous nucleotides are encoded as zeros.
        Adapted from PROcapNet code.

        Keyword arguments:
        - seq: Sequence to encode.

        Output:
        - encoding: numpy array of shape (len(seq), 4)
        '''
        ambiguous_nucs = ["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B", "X", "N"]

        sequence = seq.upper()
        if isinstance(sequence, str):
            sequence = list(sequence)

        alphabet = ["A", "C", "G", "T"]
        alphabet_lookup = {char: i for i, char in enumerate(alphabet)}

        ohe = np.zeros((len(sequence), len(alphabet)), dtype="int8")
        for i, char in enumerate(sequence):
            if char in alphabet:
                idx = alphabet_lookup[char]
                ohe[i, idx] = 1
            else:
                assert char in ambiguous_nucs, char

        return ohe

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
    