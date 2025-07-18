
import os

import pandas as pd

from Bio import SeqIO

from .GenomicElements import GenomicElements
from .BedTable import BedTable3

class ExogeneousSequences(GenomicElements):
    '''
    Class for Exogeneous sequences.

    This class inherits from GenomicElements 
    since Geneomic Elements is a special case 
    of Exogeneous sequences. Exogeneous sequences 
    class implement extra set of fasta manipulation 
    functions.

    Usually exogeneous sequences are small sets 
    of sequences compared to reference genome. 
    As a result this class read sequences into 
    genome for further analysis.
    '''
    def __init__(self, fasta_path):
        '''
        Initialize the ExogeneousSequences object.

        Key arguments:
            - fasta_path: path to the fasta file
        '''
        self.region_file_type = "bed3"
        self._anno_arr_dict = {}
        self._anno_length_dict = {}
        self.fasta_path = fasta_path

        self.sequence_df = self.read_fasta_sequences(fasta_path)
    
    @property
    def region_path(self):
        raise NotImplementedError("ExogeneousSequences does not use region_path. ")
    
    def get_region_bed_table(self):
        '''
        Return a bed table object for the sequences.
        Creates a BedTable3 with each sequence as a separate region.
        '''
        bt = BedTable3(enable_sort=False)
        
        # Create a dataframe with each sequence as a region
        # Use sequence name as chromosome, start=0, end=sequence length
        region_df = pd.DataFrame({
            "chrom": self.sequence_df.index,
            "start": 0,
            "end": [len(seq) for seq in self.sequence_df["seqs"]], 
        })
        
        bt.load_from_dataframe(region_df)
        return bt
    
    def read_fasta_sequences(self, fasta_path):
        '''
        Read fasta sequences and return a dataframe.
        '''

        seq_names = []
        seqs = []

        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_names.append(str(record.id))
            seqs.append(str(record.seq))
        
        return pd.DataFrame({"seqs": seqs,
                             }, 
                             index=seq_names,
                             )

    def get_sequence_ids(self):
        '''
        Return an array of sequence ids.
        '''
        return self.get_region_bed_table().get_chrom_names()
    
    def get_all_region_seqs(self):
        return self.sequence_df["seqs"].tolist()

    @staticmethod
    def set_parser_genome(parser):
        parser.add_argument("--fasta", 
                            help="Path to the sequence fasta file.",
                            required=True,
                            )

    @staticmethod
    def set_parser_exogeneous_sequences(parser):
        ExogeneousSequences.set_parser_genome(parser)
    
    @staticmethod
    def write_sequences_to_fasta(seq_ids: list, sequences: list, fasta_path: str):
        '''
        Write sequences to a fasta file.

        Key arguments:
            - seq_ids: list of sequence ids
            - sequences: list of sequences
            - fasta_path: path to the fasta file
        '''
        if os.path.exists(fasta_path):
            raise ValueError(f"File {fasta_path} already exists.")
        
        with open(fasta_path, "w") as f:
            for seq_id, seq in zip(seq_ids, sequences):
                f.write(f">{seq_id}\n{seq}\n")


