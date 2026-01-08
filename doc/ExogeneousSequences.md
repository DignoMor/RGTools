---
title: ExogeneousSequences Module
description: Handle non-genomic FASTA sequences with annotations and filtering
---

# ExogeneousSequences Module

Manage small sets of non-genomic sequences (FASTA), expose them as BED-like regions, and reuse `GeneralElements` helpers for annotations and per-sequence operations.

## API Reference

### ExogeneousSequences Class

**Initialization:**
```python
ExogeneousSequences(fasta_path: str)
```
- `fasta_path`: Path to the FASTA file containing the sequences.

**Key Properties:**
- `fasta_path` (str): Path to the FASTA used for sequence extraction.
- `region_file_type` (str): Always `"bed3"` for synthetic regions over sequences.
- `region_file_path` (str): Not used; accessing raises `NotImplementedError`.

**Essential Methods:**

- `get_region_bed_table() -> BedTable3`
  - Creates a mock BED3 table where each sequence is a 
    region (`chrom` = sequence id, `start` = 0, `end` = sequence length).

- `get_sequence_ids() -> list[str]`
  - Returns sequence identifiers.

- `get_all_region_seqs() -> list[str]`
  - Returns all sequences as strings.

- `get_all_region_lens() -> list[int]`
  - Returns lengths of all sequences.

- `apply_logical_filter(logical: np.ndarray, new_fasta_path: str) -> ExogeneousSequences`
  - Filters sequences by a boolean mask, writes a new FASTA, and returns a new `ExogeneousSequences`.
  - Carries over annotations filtered by the same mask.

**Inherited Helpers (see `GeneralElements.md`):**
- `get_num_regions()`
- Annotation loaders/savers (`load_region_anno_from_npy`, `load_region_anno_from_arr`, `save_anno_npy`, `save_anno_npz`)
- Annotation accessors (`get_anno_dim`, `get_anno_arr`, `get_anno_type`, `get_region_anno_by_index`, `get_region_lens`)
- `get_all_region_one_hot()` (requires length-homogeneous sequences)
- `one_hot_encoding()`

**CLI Parser Helpers:**
- `set_parser_exogeneous_sequences(parser)`: adds `--fasta`

## Examples

### Load sequences and list IDs
```python
from RGTools.ExogeneousSequences import ExogeneousSequences

es = ExogeneousSequences(fasta_path="oligos.fa")
print(es.get_sequence_ids())
print(es.get_all_region_lens())
```

### Filter sequences and keep annotations aligned
```python
import numpy as np

# Assume a stat annotation of shape (N,) or (N,1)
es.load_region_anno_from_arr("gc", np.array([0.4, 0.6, 0.55]))

mask = es.get_anno_arr("gc").reshape(-1,) > 0.5
filtered = es.apply_logical_filter(mask, "filtered_oligos.fa")

print(filtered.get_num_regions())
print(filtered.get_anno_type("gc"))  # "stat"
```

## Important Notes

- Sequences are read into memory on init; suitable for small/medium FASTA sets.
- `region_file_path` is intentionally unsupported for this class.
- Annotation shape rules (same as `GeneralElements`):
  - Stats: accept `(N,)` or `(N, 1)` and are stored as `(N, 1)`.
  - Tracks: must be `(N, L)` with `L == max(sequence length)` (no padding).
  - `get_region_anno_by_index` slices tracks to each sequence’s length; stats return scalars.
- `get_all_region_one_hot` requires length-homogeneous sequences; otherwise it raises.

