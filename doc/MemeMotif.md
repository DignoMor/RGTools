---
title: MemeMotif Module
module: MemeMotif
description: Parse, manipulate, and write MEME motif format files
version: 1.0
---

# MemeMotif Module

Parse MEME motif files, score sequences against Position Weight Matrices (PWMs), and search for motif matches in sequences.

## API Reference

### MemeMotif Class

**Initialization:**
```python
MemeMotif(file_path: str = None)
```
- `file_path`: Path to MEME file. If None, creates empty object for building new motifs.

**Key Properties:**
- `file_path` (str): Path to the original file.
- `motifs` (list): List of motif names.
- `version` (str): Version number as read from the MEME motif file.
- `alphabet` (str): Alphabet string (e.g., "ACGT").
- `strands` (list): Strand definition for the MEME motif file.
- `bg_freq` (list): Background frequencies list (same order as alphabet).
- `motif_info_dict` (dict): Dictionary of motif info. 

**Essential Methods:**

- `get_motif_list() -> list[str]`
  - Returns list of motif names

- `get_motif_pwm(motif_name: str) -> np.array`
  - Returns 2D numpy array: shape (motif_length, alphabet_length)
  - Each row is a position, each column is a base probability

- `get_motif_length(motif_name: str) -> int`
  - Returns motif length (number of positions)

- `get_alphabet() -> str`
  - Returns alphabet string

- `get_bg_freq() -> list[float]`
  - Returns background frequencies (same order as alphabet)

- `add_motif(motif_name: str, motif_info: dict) -> None`
  - Add motif. `motif_info` must contain:
    - `"alphabet_length"` (int)
    - `"motif_length"` (int)
    - `"num_source_sites"` (int)
    - `"source_eval"` (float)
    - `"pwm"` (np.array): Must be normalized (rows sum to 1.0)

- `write_meme_file(file_path: str) -> None`
  - Write MEME format file

- `clone_empty() -> MemeMotif`
  - Clone with same metadata but no motifs

**Static Methods:**

- `calculate_pwm_score(seq: str, pwm: np.array, alphabet: str = "ACGT", bg_freq: list = None, reverse_complement: bool = False) -> float`
  - Score sequence against PWM
  - **Critical**: `len(seq)` must equal `pwm.shape[0]`
  - Returns log-odds score (higher = better match)
  - If `bg_freq` is None, uses uniform frequencies

- `search_one_motif(seq: str, motif_alphabet: str, motif_pwm: np.array, bg_freq: list = None, strand: bool = False) -> np.array`
  - Search sequence for motif matches
  - Returns array of scores for each position
  - Returned Array length = sequence length
  - returned array[i] represents matching score of seq[i,i+l]
  - Positions where motif doesn't fit (at the end) are set to minimum score
  - strand is one of "+", "-" and "both". When "both" are given, 
    returned values are the higher score between fwd and rc strand.

## Important Notes

### PWM Normalization
When creating PWMs manually, rows must sum to 1.0:
```python
pwm = pwm / pwm.sum(axis=1, keepdims=True)
```

### Sequence Length Matching
`calculate_pwm_score` requires exact length match:
```python
# Check length first
if len(seq) != motif.get_motif_length(motif_name):
    # Handle mismatch
```

### Background Frequencies
- Default is uniform (0.25 for each base in ACGT)
- Order must match alphabet order
- For "ACGT", bg_freq should be [freq_A, freq_C, freq_G, freq_T]

### Score Interpretation
- Scores are log-odds (log10 ratios)
- Higher scores = better matches
- Negative scores = worse than background
- Typical good matches: > 2.0

### Search Function Behavior
- `search_one_motif` pads ends with minimum score
- For motif length L, last L-1 positions get minimum score
- Use `np.argmax(scores)` to find best position

## Dependencies

- numpy
- RGTools.utils (for reverse_complement)
