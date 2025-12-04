#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

if [ ! -d "build" ]; then
    mkdir build
fi

cd build

echo "Configuring..."
cmake -G "Unix Makefiles" ..

echo "Building..."
make

echo "Build complete: build/bin/MES_Project"
