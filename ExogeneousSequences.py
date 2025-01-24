
import pandas as pd

from Bio import SeqIO

from .GenomicElements import GenomicElements

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
    def __init__(self, region_path, region_file_type, fasta_path):
        super().__init__(region_path, region_file_type, fasta_path)

        self.sequence_df = self.read_fasta_sequences(fasta_path)
    
    
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

    @staticmethod
    def set_parser_genome(parser):
        parser.add_argument("--fasta", 
                            help="Path to the sequence fasta file.",
                            required=True,
                            )

    @staticmethod
    def set_parser_exogeneous_sequences(parser):
        ExogeneousSequences.set_parser_genome(parser)
        ExogeneousSequences.set_parser_genomic_element_region(parser)

