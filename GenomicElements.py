
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