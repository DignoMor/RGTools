"""
# Regulatory Genome Tools (RGTools)

`RGTools` is a Python package designed for efficient handling of regulatory genomic data, 
including BED-like regions, sequence extraction, motif analysis, and signal tracks.

## Overview

- **GenomicElements**: Main interface for working with genomic regions and reference genomes.
- **GeneralElements**: Base classes and shared logic for element collections.
- **ExogeneousSequences**: Specialized handling for sequences outside the reference genome.
- **MemeMotif**: Parser and scorer for MEME-formatted motifs.
- **BedTable**: Utilities for loading and manipulating BED files as pandas DataFrames.
- **ListFile**: Simple utility for reading and handling single-column list files.
"""

# Export modules to ensure they appear in the documentation submodules list
from . import GenomicElements as GenomicElements_mod
from . import ExogeneousSequences as ExogeneousSequences_mod
from . import MemeMotif as MemeMotif_mod
from . import BedTable as BedTable_mod
from . import ListFile as ListFile_mod
from . import BwTrack as BwTrack_mod
from . import utils as utils_mod
from . import SNP_utils as SNP_utils_mod
from . import GTF_utils as GTF_utils_mod
from . import exceptions as exceptions_mod
from . import logging as logging_mod
from . import GeneralElements as GeneralElements_mod

# Export main classes and modules for easier access
from .GenomicElements import GenomicElements
from .ExogeneousSequences import ExogeneousSequences
from .MemeMotif import MemeMotif
from .BedTable import BedTable3, BedTable6, BedTable6Plus, BedTable3Plus, BedTablePairEnd
from .ListFile import ListFile
from .BwTrack import SingleBwTrack, PairedBwTrack

__all__ = [
    # Submodules
    "GenomicElements",
    "ExogeneousSequences",
    "MemeMotif",
    "BedTable",
    "ListFile",
    "BwTrack",
    "utils",
    "SNP_utils",
    "GTF_utils",
    "exceptions",
    "logging",
    "GeneralElements",
    # Top-level classes for convenience
    "BedTable3",
    "BedTable6",
    "BedTable6Plus",
    "BedTable3Plus",
    "BedTablePairEnd",
    "SingleBwTrack",
    "PairedBwTrack"
]

