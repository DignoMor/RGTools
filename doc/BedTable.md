---
title: BedTable Module
description: Load and manipulate BED files as pandas DataFrames
---

# BedTable Module

The `BedTable` module provides classes for loading, validating, and manipulating various BED file formats. It uses `pandas.DataFrame` for efficient data storage and manipulation, while providing a genomic-aware interface.

## API Reference

### Core Classes

The module provides specialized classes for different BED flavors:

- **`BedTable3`**: Standard 3-column BED (`chrom`, `start`, `end`).
- **`BedTable6`**: Standard 6-column BED (adds `name`, `score`, `strand`).
- **`BedTable3Plus` / `BedTable6Plus`**: Base classes for BED files with additional custom columns.
- **`BedTablePairEnd`**: Support for paired-end genomic regions (`chrom`, `start`, `end`, `chrom2`, `start2`, `end2`, etc.).

### Initialization

Most classes accept an `enable_sort` parameter (default `True`). If enabled, the table will be automatically sorted by chromosome and position upon loading.

```python
from RGTools.BedTable import BedTable3
bt = BedTable3(enable_sort=True)
```

### Key Methods

#### Data I/O
- **`load_from_file(ipath)`**: Load data from a BED file. Supports `stdin`.
- **`load_from_dataframe(df, column_map=None)`**: Load data from an existing pandas DataFrame.
- **`write(opath)`**: Save data to a BED file. Supports `stdout`.
- **`to_dataframe()`**: Returns a copy of the underlying `pandas.DataFrame`.

#### Filtering and Subsetting
- **`apply_logical_filter(logical_array)`**: Returns a new table containing only the regions where the boolean array is `True`.
- **`region_subset(chrom, start, end)`**: Returns a new table with regions fully contained within the specified coordinates.
- **`subset_by_index(index_array)`**: Returns a new table containing regions at the specified integer indices.

#### Searching and Inspection
- **`search_region(chrom, start, end, overlapping_base=1)`**: Returns indices of regions that overlap with the query.
- **`get_region_by_index(index)`**: Returns a `BedRegion` object for the specified index.
- **`iter_regions()`**: Returns an iterator over regions (returns `BedRegion` objects).
- **`__len__`**: Use `len(bt)` to get the number of regions.

#### Accessors (Numpy Arrays)
- **`get_chrom_names()`**, **`get_start_locs()`**, **`get_end_locs()`**: Access basic coordinates as numpy arrays.
- **`get_region_names()`**, **`get_region_scores()`**, **`get_region_strands()`**: Access basic bed6 information.
- **`get_pair_names()`**, **`get_pair_scores()`**: Access paired-end specific information (for `BedTablePairEnd`).
- **`get_region_extra_column(column_name)`**: (For `Plus` classes) Access custom columns as numpy arrays.

### Helper Factory Methods (in `GenomicElements`)

While not separate classes, common BED formats are supported via factory methods in `GenomicElements`:
- `BedTableNarrowPeak()`: For ENCODE narrowPeak format.
- `BedTableTREBed()`: For TRE-specific BED formats.
- `BedTableBedGraph()`: For bedGraph files.

## Example

```python
from RGTools.BedTable import BedTable6

# Load and search
bt = BedTable6()
bt.load_from_file("peaks.bed")

# Find regions overlapping a specific locus
indices = bt.search_region("chr1", 1000, 2000)
for idx in indices:
    region = bt.get_region_by_index(idx)
    print(f"Found {region.name} at {region.start}-{region.end}")

# Export filtered table
high_score = bt.apply_logical_filter(bt.to_dataframe()["score"] > 500)
high_score.write("top_peaks.bed")
```
