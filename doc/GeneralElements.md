---
title: GeneralElements Module
description: Shared base for genomic/exogenous element collections and annotations
---

# GeneralElements Module

Abstract base class for `GenomicElements` and `ExogeneousSequences`. Provides shared sequence retrieval, one-hot encoding, annotation management, and filtering hooks that subclasses extend with their own region/sequence implementations.

## API Reference

### GeneralElements Class

**Initialization:**

- Abstract base; instantiate through a subclass (e.g., `GenomicElements`).

**Abstract Properties:**

- `fasta_path` (str): Path to the FASTA used for sequence extraction.
- `region_file_type` (str): Region file type string used by the subclass.
- `region_file_path` (str): Path to the region file.

**Abstract Methods:**
- `get_region_bed_table() -> BedTable`
  - Returns the underlying bed table object for the subclass.
- `apply_logical_filter(logical: np.ndarray, new_path: str) -> GeneralElements`
  - Filters regions using a boolean mask and writes the filtered regions to `new_path`.

**Essential Methods:**

- `get_num_regions() -> int`
  - Returns the number of regions (delegates to the bed table length).

- `get_all_region_seqs() -> list[str]`
  - Loads the entire FASTA into memory and returns sequences for all regions.
  - Raises if a chromosome in the bed table is missing from the FASTA.

- `get_region_lens() -> np.ndarray`
  - Returns an array of per-region lengths (`end - start`).

- `get_anno_type(anno_name: str) -> str`: 
  - Return the type of the annotation (track or stat)

- `load_region_anno_from_npy(anno_name: str, npy_path: str) -> None`
  - Loads per-region annotation from `.npy` **or** single-array `.npz`.

- `load_region_anno_from_arr(anno_name: str, anno_arr: np.ndarray) -> None`
  - Loads per-region annotation from a numpy array.
  - `anno_arr.shape[0]` must equal `get_num_regions()`.
  - Stats accept `(N,)` or `(N, 1)` and are stored as `(N, 1)`.
  - Track annotations accept only `(N, L), L = max(region length)`. A Value error will be 
    raised otherwise.

- `load_region_track_from_list(anno_name: str, anno_list: list[np.ndarray]) -> None`
  - Loads per-region annotation track from a python list.
  - useful when loading annotations for non-length-homogeneous elements
  - `anno_list[i]` must have the same length as the ith element

- `get_anno_arr(anno_name: str) -> np.ndarray`
  - Returns the stored annotation array.

- `get_region_lens() -> np.ndarray`
  - Return the lengths of each region as an array.

- `get_region_anno_by_index(anno_name: str, index: int)`
  - Returns annotation for a specific region.
  - For padded annotations, slices to the region’s length.

- `save_anno_npy(anno_name: str, npy_path: str) -> None`
  - Saves annotation to `.npy`.

- `save_anno_npz(anno_name: str, npz_path: str) -> None`
  - Saves annotation to compressed `.npz`.

**Methods restricted to length-homogeneous Elements**

- `get_all_region_one_hot() -> np.ndarray`
  - Returns array of shape `(num_regions, region_length, 4)`.
  - Requires **length-homogeneous regions**; otherwise raises `ValueError` listing lengths.

**Static Methods**

- `one_hot_encoding(seq: str) -> np.ndarray`
  - Static method. Returns `(len(seq), 4)` int8 array for A/C/G/T.
  - Ambiguous nucleotides (e.g., `N`, `R`, `Y`, etc.) are encoded as all zeros.

## Examples

### Use via a subclass (GenomicElements)
```python
from RGTools.GenomicElements import GenomicElements

ge = GenomicElements(
    region_file_path="peaks.bed3",
    region_file_type="bed3",
    fasta_path="hg38.fa",
)

# Base-class helpers
seqs = ge.get_all_region_seqs()
one_hot = ge.get_all_region_one_hot()  # shape: (N, L, 4)

import numpy as np
anno = np.random.randn(ge.get_num_regions(), 3)
ge.load_region_anno_from_arr("my_anno", anno)
ge.save_anno_npz("my_anno", "my_anno.npz")

mask = ge.get_anno_arr("my_anno")[:, 0] > 0
filtered = ge.apply_logical_filter(mask, "filtered.bed3")
print(filtered.get_num_regions())
```

## Important Notes

### Lengh-homogeneous and non-length-homogeneous regions

A set of Elements are lengh-homogeneous if everyone of them 
are of the same length. Otherwise it is non-length-homogeneous.

- 2D annotations are padded with zeros to the longest region length.
- `get_anno_arr` returns the padded array; `get_region_anno_by_index` slices per region length.
- `get_all_region_one_hot` still requires length-homogeneous regions.

### Annotation types

There are 2 types of annotations

- track: store signal tracks for each element. The second dimension is the element length.
  - for length-homogeneous elements, the second dimension is the length.
  - for non-length-homogeneous elements, the second dimension varies by elements. 
    However, to store the track as np array we pad them to match the longest 
    element and slice the track at the time of retrival.
- stat: store statistics for each element, the second dimension is of size 1

### Other Notes 

- `get_all_region_seqs()` loads the FASTA into memory; use when the genome fits in memory.
- `get_all_region_one_hot()` additionally materializes a `(N, L, 4)` array; ensure consistent region lengths first.
- `.npz` loading via `load_region_anno_from_npy()` only works when the file has exactly one array; otherwise a `ValueError` lists available keys.
- One-hot encoding treats ambiguous nucleotides as zeros.
  



