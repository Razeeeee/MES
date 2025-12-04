#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

if [ ! -f "build/bin/MES_Project" ]; then
    echo "Executable not found. Run 'make build' first."
    exit 1
fi

./build/bin/MES_Project
