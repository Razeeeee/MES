#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root directory
cd "$PROJECT_ROOT"

# Build the project first
./scripts/build.sh

# Check if build was successful
if [ $? -eq 0 ]; then
    echo ""
    ./build/bin/MES_Project
else
    echo "Build failed."
    exit 1
fi
