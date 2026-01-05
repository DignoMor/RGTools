#!/bin/bash
# Generate documentation using pdoc

set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the project root (parent of scripts/)
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

cd "$PROJECT_ROOT"

# Check if pdoc is installed
if ! command -v pdoc &> /dev/null; then
    echo "pdoc is not installed. Installing..."
    pip install pdoc
fi

echo "Generating documentation..."
pdoc -o github_pages $PROJECT_ROOT

echo "Documentation generated in github_pages/ folder"

