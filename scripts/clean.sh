#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root directory
cd "$PROJECT_ROOT"

# Remove build directory and its contents
if [ -d "build" ]; then
    rm -rf build
    echo "Cleaned build directory."
else
    echo "Nothing to clean."
fi
