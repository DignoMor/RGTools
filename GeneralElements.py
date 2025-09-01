import abc
import numpy as np

from Bio import SeqIO

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
