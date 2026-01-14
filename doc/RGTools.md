---
title: RGTools Package
description: Regulatory Genome Tools - A collection of utilities for genomic data analysis
---

# Regulatory Genome Tools (RGTools)

`RGTools` is a Python package designed for efficient handling of regulatory genomic data, including BED-like regions, sequence extraction, motif analysis, and signal tracks.

## Overview

The package is organized into several key modules:

- **[GenomicElements](GenomicElements.md)**: Main interface for working with genomic regions and reference genomes.
- **[GeneralElements](GeneralElements.md)**: Base classes and shared logic for element collections.
- **[ExogeneousSequences](ExogeneousSequences.md)**: Specialized handling for sequences outside the reference genome.
- **[MemeMotif](MemeMotif.md)**: Parser and scorer for MEME-formatted motifs.
- **[BedTable](BedTable.md)**: Utilities for loading and manipulating BED files as pandas DataFrames.
- **[ListFile](ListFile.md)**: Simple utility for reading and handling single-column list files.

## Key Features

- Efficient one-hot encoding of genomic sequences.
- Unified interface for various BED formats (bed3, bed6, narrowPeak, etc.).
- Robust annotation management (tracks and statistics) linked to genomic regions.
- Fast sequence retrieval from large FASTA files.
