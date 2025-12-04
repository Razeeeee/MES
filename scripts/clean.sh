#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

if [ -d "build" ]; then
    rm -rf build
    echo "Cleaned build directory."
else
    echo "Nothing to clean."
fi
