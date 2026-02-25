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
  - Filters regions using a **mask** and writes the filtered regions to `new_path`.

**Essential Getter Methods:**

- `get_num_regions() -> int`
  - Returns the number of regions (delegates to the bed table length).

- `get_all_region_seqs() -> list[str]`
  - Loads the entire FASTA into memory and returns sequences for all regions.
  - Raises if a chromosome in the bed table is missing from the FASTA.

- `get_region_lens() -> np.ndarray`
  - Returns an array of per-region lengths (`end - start`).

- `get_anno_type(anno_name: str) -> str`: 
  - Return the type of the annotation (`track`, `stat`, `mask`, or `array`).

- `get_anno_dim(anno_name: str) -> int | tuple[int, ...]`
  - Return stored annotation dimension metadata.
  - `track`: `max_region_len`; `stat`/`mask`: `1`; `array`: trailing shape tuple.

The following methods allow retrieval of all regions' annotations.

- `get_track_list(anno_name: str) -> list[np.ndarray]`
  - Return track annotations as a list.
  - Each returned track has shape `(region_length,)`.

- `get_stat_arr(anno_name: str) -> np.ndarray`
  - Return stat annotation array with shape `(N, 1)`.

- `get_mask_arr(anno_name: str) -> np.ndarray`
  - Return mask annotation array with shape `(N, 1)`.

- `get_arr_anno(anno_name: str) -> np.ndarray`
  - Return array annotation with shape `(N, ...)`.

The following methods allow index-based retrieval of 
one region's annotation.

- `get_region_track_by_index(anno_name: str, index: int) -> np.ndarray`
  - Return track annotation for a specific region.
  - Return shape: `(region_size_index,)`.

- `get_region_stat_by_index(anno_name: str, index: int) -> numerical`
  - Return stat annotation for a specific region.

- `get_region_mask_by_index(anno_name: str, index: int) -> bool`
  - Return mask annotation for a specific region.

- `get_region_array_by_index(anno_name: str, index: int) -> np.ndarray`
  - Return array annotation for a specific region.

**Essential Setter Methods:**

- `load_region_track_from_list(anno_name: str, anno_list: Iterable[np.ndarray]) -> None`
  - Loads per-region annotation **track** from a python list.
  - `anno_list[i]` must have the same length as the ith element

- `load_region_stat_from_arr(anno_name: str, anno_arr: np.ndarray) -> None`
  - Loads per region **stat** from a numpy array.

- `load_mask_from_arr(anno_name:str, anno_arr: np.ndarray) -> None`
  - Loads a region **mask** from a numpy array.

- `load_region_array_from_arr(anno_name: str, anno_arr: np.ndarray) -> None`
  - Load region **array** annotations.
  - `anno_arr` must have shape `(n_region, anno_shape...)`.
  - For example, one-hot array annotation has shape `(n_region, region_length, 4)`.

**Essential IO Methods:**

- `save_anno_npy(anno_name: str, npy_path: str) -> None`
  - Saves annotation to `.npy`.

- `save_anno_npz(anno_name: str, npz_path: str) -> None`
  - Saves annotation to compressed `.npz`.

- `load_region_anno_from_npy(anno_name: str, npy_path: str, anno_type: str = "array") -> None`
  - Loads per-region annotation from `.npy` **or** single-array `.npz`.
  - When loading tracks, set `anno_type="track"` so track semantics are preserved.

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
stat = np.random.randn(ge.get_num_regions())
ge.load_region_stat_from_arr("my_stat", stat)
ge.save_anno_npz("my_stat", "my_stat.npz")

mask = ge.get_stat_arr("my_stat")[:, 0] > 0
filtered = ge.apply_logical_filter(mask, "filtered.bed3")
print(filtered.get_num_regions())
```

## Important Notes

### Length-homogeneous and non-length-homogeneous regions

A set of Elements are length-homogeneous if every one of them 
are of the same length. Otherwise it is non-length-homogeneous.

- 2D track annotations are padded with zeros to the longest region length.
- `get_track_list` and `get_region_track_by_index` slice track values to each region length.
- `get_all_region_one_hot` still requires length-homogeneous regions.

### Annotation types

There are 4 annotation types:

- **track**: signal track per element.
  - Setter: `load_region_track_from_list(anno_name, anno_list)`
  - Getter(s): `get_track_list(anno_name)`, `get_region_track_by_index(anno_name, index)`
  - Stored as shape `(N, max_region_len)` with zero-padding for non-length-homogeneous regions.
- **stat**: scalar statistic per element.
  - Setter: `load_region_stat_from_arr(anno_name, anno_arr)`
  - Getter(s): `get_stat_arr(anno_name)`, `get_region_stat_by_index(anno_name, index)`
  - Stored as shape `(N, 1)`.
- **mask**: boolean scalar per element.
  - Setter: `load_mask_from_arr(anno_name, anno_arr)`
  - Getter(s): `get_mask_arr(anno_name)`, `get_region_mask_by_index(anno_name, index)`
  - Stored as shape `(N, 1)`.
- **array**: fixed-shape array payload per element.
  - Setter: `load_region_array_from_arr(anno_name, anno_arr)`
  - Getter(s): `get_arr_anno(anno_name)`, `get_region_array_by_index(anno_name, index)`
  - Stored as shape `(N, ...)` where all elements share the same trailing shape.


Interally, all annotations are stored as numpy array. To avoid 
array dimension issues with non-length-homogeneous elements, 
tracks are padded with zeros to the maximum region size before 
storing.

### Other Notes 

- `get_all_region_seqs()` loads the FASTA into memory; use when the genome fits in memory.
- `get_all_region_one_hot()` additionally materializes a `(N, L, 4)` array; ensure consistent region lengths first.
- `.npz` loading via `load_region_anno_from_npy()` only works when the file has exactly one array; otherwise a `ValueError` lists available keys.
- One-hot encoding treats ambiguous nucleotides as zeros.
- For general usage, region tracks are loaded from python lists for 
  the best compatibility. Masks and stats are loaded from arr 
  since the data can be stored in ndarray of known shape.



