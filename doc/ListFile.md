---
title: ListFile Module
description: Read and handle simple single-column list files
---

# ListFile Module

The `ListFile` module provides a simple utility for reading and handling text files that contain one item per line, such as lists of gene symbols, region IDs, or sample names.

## API Reference

### ListFile Class

**Initialization:**
```python
ListFile(filter_empty_lines=True)
```
- `filter_empty_lines`: If `True`, whitespace-only lines are ignored.

**Methods:**

- **`read_file(file_path)`**: Load items from a file. Supports `stdin`.
- **`get_contents(dtype="str")`**: Returns the items as a numpy array. 
- **`get_num_lines()`**: Returns the number of items loaded.

## Example

```python
from RGTools.ListFile import ListFile

lf = ListFile()
lf.read_file("genes.txt")

print(f"Loaded {lf.get_num_lines()} genes.")
genes = lf.get_contents()
for gene in genes:
    print(gene)
```
