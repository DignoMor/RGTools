---
title: GenomicElements Module
description: Work with genomic regions (BED-like files) and a reference genome FASTA
---

# GenomicElements Module

Load genomic regions from BED-like files, retrieve reference sequences from a genome FASTA, create one-hot encodings, attach per-region annotations, and filter/export regions.

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
- `region_file_type` (str): Region file type string (e.g., `"bed3"`).
- `fasta_path` (str): Path to the genome FASTA used for sequence extraction.

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
  - Returns a **copy** of the underlying bed table object (e.g., `BedTable3`, `BedTable6Plus`).

- `get_region_seq(chrom: str, start: int, end: int) -> str | None`
  - Returns the reference sequence for the given region (BED coordinates).
  - Returns `None` if `chrom` is not found in the FASTA.

- `export_exogeneous_sequences(fasta_path: str) -> None`
  - Writes each region as a FASTA record to `fasta_path`.
  - Raises if `fasta_path` already exists.

- `apply_logical_filter(logical: np.ndarray, new_region_file_path: str) -> GenomicElements`
  - Filters regions using a boolean mask and writes filtered regions to `new_region_file_path`.
  - Also filters and carries over any loaded annotations.

**Inherited Helpers (see `GeneralElements.md`):**
- `get_num_regions()`
- `get_all_region_seqs()`
- `get_all_region_one_hot()`
- Annotation loaders/savers (`load_region_anno_from_npy`, `load_region_anno_from_arr`, `save_anno_npy`, `save_anno_npz`)
- Annotation accessors (`get_anno_dim`, `get_anno_arr`, `get_anno_type`, `get_region_anno_by_index`, `get_region_lens`)
- `one_hot_encoding()`

**CLI Parser Helpers:**
- `set_parser_genome(parser)`: adds `--fasta_path`
- `set_parser_genomic_element_region(parser)`: adds `--region_file_path` and `--region_file_type`

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
mask = ge.get_anno_arr("my_anno")[:, 0] > 0
filtered = ge.apply_logical_filter(mask, "filtered.bed3")
print(filtered.get_num_regions())
```

## Important Notes

### Coordinate Convention
All region coordinates use **BED convention**: 0-based, end-exclusive slicing (`seq[start:end]`).

### FASTA Chromosome IDs Must Match
`chrom` values in your region file must match FASTA record IDs (e.g., `chr1` vs `1`).

### Memory Considerations
- `get_region_seq()` scans the FASTA until it finds the matching record.

### Annotation Shape Rules (see `GeneralElements.md`)
- Stats: accept `(N,)` or `(N, 1)` and are stored as `(N, 1)`.
- Tracks: must be `(N, L)` with `L == max(region length)` (no padding performed).
- `get_region_anno_by_index` slices tracks to each region’s length; stats return scalars.
- `get_all_region_one_hot` requires length-homogeneous regions; otherwise it raises.

For full annotation, one-hot encoding, and `.npz` handling details, see `GeneralElements.md`.