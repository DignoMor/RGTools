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
- `region_file_path`: Path to the region file (BED-like).
- `region_file_type`: One of the supported region formats (see below).
- `fasta_path`: Path to a genome FASTA file. FASTA record IDs must match the `chrom` field in the region file.

**Key Properties:**
- `fasta_path` (str): Path to the genome FASTA file.
- `region_file_type` (str): Type of the region file.
- `region_file_path` (str): Path to the original region file.

**Supported `region_file_type` values:**
`GenomicElements.get_region_file_suffix2class_dict()` supports:
- `bed3`: 3-column BED (`chrom`, `start`, `end`)
- `bed6`: 6-column BED (`chrom`, `start`, `end`, `name`, `score`, `strand`)
- `bed3gene`: BED3 + `gene_symbol`
- `bed6gene`: BED6 + `gene_symbol`
- `narrowPeak`: BED6 + (`signalValue`, `pValue`, `qValue`, `peak`)
- `TREbed`: BED3 + (`name`, `fwdTSS`, `revTSS`)
- `bedGraph`: BED3 + (`dataValue`)

**Essential Methods:**

- `get_region_bed_table() -> BedTable`
  - Returns a copy of the underlying bed table object (e.g., `BedTable3`, `BedTable6Plus`).
  - Can be used to iterate over regions or perform bed table operations.

- `get_region_seq(chrom: str, start: int, end: int) -> str | None`
  - Returns the reference sequence for the given region (BED coordinates).
  - Coordinates follow BED convention (0-based, half-open).
  - Returns `None` if `chrom` is not found in the FASTA.

- `export_exogeneous_sequences(fasta_path: str) -> None`
  - Writes each region as a FASTA record to `fasta_path`.
  - Raises if `fasta_path` already exists.
  - FASTA header format: `>chrom:start-end`.

- `apply_logical_filter(logical: np.ndarray, new_region_file_path: str) -> GenomicElements`
  - Filters regions using a boolean mask and writes filtered regions to `new_region_file_path`.
  - Also filters and carries over loaded annotations.

**Inherited Helpers (see `GeneralElements.md`):**
- `get_num_regions()`
- `get_all_region_seqs()`
- `get_all_region_one_hot()`
- Annotation loaders/savers (`load_region_anno_from_npy`, `load_region_track_from_list`, `load_region_stat_from_arr`, `load_mask_from_arr`, `load_region_array_from_arr`, `save_anno_npy`, `save_anno_npz`)
- Annotation accessors (`get_anno_type`, `get_track_list`, `get_stat_arr`, `get_mask_arr`, `get_arr_anno`, `get_region_track_by_index`, `get_region_stat_by_index`, `get_region_mask_by_index`, `get_region_array_by_index`, `get_region_lens`)
- `one_hot_encoding()`

**CLI Parser Helpers:**
- `set_parser_genome(parser)`: adds `--fasta_path`
- `set_parser_genomic_element_region(parser)`: adds `--region_file_path` and `--region_file_type`

**Additional Static Helpers:**
- `BedTable6Gene(enable_sort: bool = True) -> BedTable6Plus`
  - Build bed table parser for `bed6gene` with extra `gene_symbol` field.
- `BedTable3Gene(enable_sort: bool = True) -> BedTable3Plus`
  - Build bed table parser for `bed3gene` with extra `gene_symbol` field.
- `BedTableNarrowPeak(enable_sort: bool = True) -> BedTable6Plus`
  - Build parser for `narrowPeak` with `signalValue`, `pValue`, `qValue`, `peak`.
- `BedTableBedGraph(enable_sort: bool = True) -> BedTable3Plus`
  - Build parser for `bedGraph` with extra `dataValue` field.
- `BedTableTREBed(enable_sort: bool = True) -> BedTable3Plus`
  - Build parser for `TREbed` with `name`, `fwdTSS`, `revTSS`.
- `merge_genomic_elements(left_ge, right_ge, output_region_path, anno2merge, sort_new_ge=True) -> GenomicElements`
  - Merge two `GenomicElements`, merge selected annotations in output order, and write merged regions to `output_region_path`.

## Examples

### Initialize and access regions
```python
from RGTools.GenomicElements import GenomicElements

ge = GenomicElements(
    region_file_path="peaks.bed3",
    region_file_type="bed3",
    fasta_path="hg38.fa",
)

bt = ge.get_region_bed_table()
print(ge.get_num_regions())  # inherited from GeneralElements
```

### Get a single region sequence
```python
seq = ge.get_region_seq("chr1", 100_000, 100_050)  # BED coords: 0-based, end-exclusive
if seq is None:
    raise ValueError("chrom not found in FASTA")
print(seq)
```

### Export regions as a FASTA
```python
ge.export_exogeneous_sequences("regions.fa")
```

### Filter regions and keep annotations aligned
```python
mask = ge.get_stat_arr("my_stat")[:, 0] > 0
filtered = ge.apply_logical_filter(mask, "filtered.bed3")
print(filtered.get_num_regions())
```

## Important Notes

### Coordinate Convention
- All coordinates use BED convention: 0-based, half-open intervals.
- `start` is inclusive and `end` is exclusive.
- Example: `chr1:100-200` includes positions 100-199.

### FASTA Chromosome IDs Must Match
`chrom` values in your region file must match FASTA record IDs (e.g., `chr1` vs `1`).

### Memory Considerations
- `get_all_region_seqs()` and `get_all_region_one_hot()` load the whole genome into memory.
- `get_region_seq()` scans the FASTA until it finds the matching record.
- For large genomes or many regions, use iterative access when possible.

### Annotation Shape Rules (see `GeneralElements.md`)
- Stats: stored as `(N, 1)`, retrieved by `get_stat_arr` / `get_region_stat_by_index`.
- Masks: stored as `(N, 1)`, retrieved by `get_mask_arr` / `get_region_mask_by_index`.
- Tracks: stored with width `max(region length)` and retrieved by `get_track_list` / `get_region_track_by_index`.
- Arrays: stored as `(N, ...)`, retrieved by `get_arr_anno` / `get_region_array_by_index`.
- `get_all_region_one_hot` requires length-homogeneous regions; otherwise it raises.

### Filtering Regions
- `apply_logical_filter()` returns a new `GenomicElements` object.
- Original object is unchanged.
- All loaded annotations are filtered to match the remaining regions.

For full annotation, one-hot encoding, and `.npz` handling details, see `GeneralElements.md`.

## Dependencies

- numpy
- BioPython (SeqIO)
- RGTools.BedTable (BedTable3, BedTable6, BedTable6Plus, BedTable3Plus)
- RGTools.GeneralElements (base class)