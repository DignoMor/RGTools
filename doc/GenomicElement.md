---
title: GenomicElements Module
description: Load, manipulate, and extract sequences from genomic region files
---

# GenomicElements Module

Load genomic region files (BED, narrowPeak, bedGraph, etc.), extract sequences from genome FASTA files, and manage annotations for genomic elements.

## API Reference

### GenomicElements Class

**Initialization:**
```python
GenomicElements(region_file_path: str, region_file_type: str, fasta_path: str)
```
- `region_file_path`: Path to the region file (BED, narrowPeak, etc.)
- `region_file_type`: Type of region file (see `get_region_file_suffix2class_dict()`)
- `fasta_path`: Path to the genome FASTA file

**Key Properties:**
- `fasta_path` (str): Path to the genome FASTA file
- `region_file_type` (str): Type of the region file
- `region_file_path` (str): Path to the original region file

**Essential Methods:**

- `get_region_bed_table() -> BedTable`
  - Returns a copy of the bed table object for the region file
  - Can be used to iterate over regions or perform bed table operations

- `get_region_seq(chrom: str, start: int, end: int) -> str`
  - Extract sequence for a specific genomic region
  - Coordinates follow BED convention (0-based, half-open)
  - Returns sequence string, or None if chromosome not found

- `get_all_region_seqs() -> list[str]`
  - Get sequences for all regions
  - **Note**: Reads entire genome file into memory
  - Returns list of sequence strings, one per region

- `get_all_region_one_hot() -> np.array`
  - Get one-hot encoding for all regions
  - **Note**: Memory intensive, reads entire genome into memory
  - Returns numpy array of shape (num_regions, sequence_length, 4)
  - All regions must have the same length
  - Ambiguous nucleotides encoded as zeros

- `export_exogeneous_sequences(fasta_path: str) -> None`
  - Export all regions as sequences in FASTA format
  - Raises ValueError if output file already exists
  - FASTA headers format: `>chrom:start-end`

- `apply_logical_filter(logical: np.array, new_region_file_path: str) -> GenomicElements`
  - Filter regions using a boolean array
  - `logical`: Boolean array of length equal to number of regions
  - Creates new GenomicElements object with filtered regions
  - Preserves annotations (filters them accordingly)

- `get_num_regions() -> int`
  - Returns the number of regions

- `load_region_anno_from_arr(anno_name: str, anno_arr: np.array) -> None`
  - Load annotation array for regions
  - `anno_arr` first dimension must match number of regions
  - Supports 1D or 2D annotation arrays

- `load_region_anno_from_npy(anno_name: str, npy_path: str) -> None`
  - Load annotation from numpy file
  - Convenience wrapper around `load_region_anno_from_arr`

- `get_anno_arr(anno_name: str) -> np.array`
  - Get annotation array by name
  - Returns the stored annotation array

- `get_anno_dim(anno_name: str) -> int`
  - Get annotation dimension (1 for 1D, second dimension size for 2D)

- `save_anno_npy(anno_name: str, npy_path: str) -> None`
  - Save annotation to numpy file

- `save_anno_npz(anno_name: str, npz_path: str) -> None`
  - Save annotation to compressed numpy file

**Static Methods:**

- `get_region_file_suffix2class_dict() -> dict`
  - Returns dictionary mapping region file types to BedTable classes
  - Valid types: `"bed3"`, `"bed6"`, `"bed6gene"`, `"bed3gene"`, `"narrowPeak"`, `"TREbed"`, `"bedGraph"`

- `BedTable6Gene(enable_sort: bool = True) -> BedTable6Plus`
  - Helper to create BedTable for bed6gene format
  - Includes extra column: `gene_symbol` (str)

- `BedTable3Gene(enable_sort: bool = True) -> BedTable3Plus`
  - Helper to create BedTable for bed3gene format
  - Includes extra column: `gene_symbol` (str)

- `BedTableNarrowPeak(enable_sort: bool = True) -> BedTable6Plus`
  - Helper to create BedTable for narrowPeak format
  - Includes extra columns: `signalValue` (float), `pValue` (float), `qValue` (float), `peak` (int)

- `BedTableBedGraph(enable_sort: bool = True) -> BedTable3Plus`
  - Helper to create BedTable for bedGraph format
  - Includes extra column: `dataValue` (float)

- `BedTableTREBed(enable_sort: bool = True) -> BedTable3Plus`
  - Helper to create BedTable for TREbed format
  - Includes extra columns: `name` (str), `fwsTSS` (int), `revTSS` (int)

- `set_parser_genome(parser) -> None`
  - Add `--fasta_path` argument to argparse parser

- `set_parser_genomic_element_region(parser) -> None`
  - Add `--region_file_path` and `--region_file_type` arguments to argparse parser

- `one_hot_encoding(seq: str) -> np.array`
  - Convert DNA sequence to one-hot encoding
  - Returns array of shape (len(seq), 4) with dtype int8
  - Ambiguous nucleotides (Y, R, W, S, K, M, D, V, H, B, X, N) encoded as zeros
  - Order: A, C, G, T

## Important Notes

### Coordinate Convention
- All coordinates follow BED convention: 0-based, half-open intervals
- `start` is inclusive, `end` is exclusive
- Example: `chr1:100-200` includes positions 100-199

### Memory Considerations
- `get_all_region_seqs()` and `get_all_region_one_hot()` load the entire genome into memory
- For large genomes or many regions, consider using `get_region_seq()` iteratively
- `get_all_region_one_hot()` requires all regions to have the same length

### Region File Types
Supported region file types:
- `bed3`: Standard 3-column BED format
- `bed6`: Standard 6-column BED format
- `bed6gene`: BED6 with gene symbol annotation
- `bed3gene`: BED3 with gene symbol annotation
- `narrowPeak`: ENCODE narrowPeak format (includes signal, p-value, q-value, peak)
- `bedGraph`: UCSC bedGraph format (includes data value)
- `TREbed`: TRE format (includes name, forward TSS, reverse TSS)

### Annotation Management
- Annotations are stored by name and can be 1D or 2D arrays
- First dimension must always match the number of regions
- Annotations are preserved when filtering regions
- Use `get_anno_dim()` to check annotation dimensionality

### Sequence Extraction
- `get_region_seq()` searches for chromosome by ID matching FASTA record IDs
- Returns `None` if chromosome not found (does not raise exception)
- For better performance with multiple queries, consider loading genome once

### Filtering Regions
- `apply_logical_filter()` creates a new GenomicElements object
- Original object is unchanged
- All annotations are automatically filtered to match remaining regions
- New region file is written to disk

## Dependencies

- numpy
- BioPython (SeqIO)
- RGTools.BedTable (BedTable3, BedTable6, BedTable6Plus, BedTable3Plus)
- RGTools.GeneralElements (base class)

