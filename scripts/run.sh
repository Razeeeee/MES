#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root directory
cd "$PROJECT_ROOT"

# Check if the executable exists
if [ ! -f "build/bin/MES_Project" ]; then
    echo "Executable not found. Run 'make build' first."
    exit 1
fi

./build/bin/MES_Project
