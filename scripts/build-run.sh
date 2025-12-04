#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

./scripts/build.sh

if [ $? -eq 0 ]; then
    echo ""
    ./build/bin/MES_Project
else
    echo "Build failed."
    exit 1
fi
