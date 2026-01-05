
# Regulatory Genome Tools (RGTools)

## Testing

Run unittest-based tests under `tests/` folder.
Before running tests, put large example files 
under `large_files/`. Files to be downloaded can 
be found in `download_large_files.sh`.

## Documentation

Documentation is automatically generated using [pdoc](https://pdoc.dev/) and hosted on GitHub Pages.

### Generating Documentation Locally

To generate documentation locally:

```bash
./scripts/generate_docs.sh
```

The documentation will be generated in the `github_pages/` folder. 

### Automatic Generation and Deployment

Documentation is automatically regenerated and deployed to GitHub Pages on every push to the `main` branch via GitHub Actions. 

