.PHONY: build run build-run clean init help

# Default target
all: build

# Build the project
build:
	@chmod +x scripts/build.sh
	@./scripts/build.sh

# Run the project (already built)
run:
	@chmod +x scripts/run.sh
	@./scripts/run.sh

# Build and run the project
build-run:
	@chmod +x scripts/build-run.sh
	@./scripts/build-run.sh

# Clean build artifacts
clean:
	@chmod +x scripts/clean.sh
	@./scripts/clean.sh

# Initialize project permissions (run this on a new machine)
init:
	@echo "Setting script permissions..."
	@chmod +x scripts/*.sh
	@echo "Done. Ready to build."

# Display help information
help:
	@echo "Available commands:"
	@echo "  make init      - Set script permissions"
	@echo "  make build     - Build project"
	@echo "  make run       - Run executable"
	@echo "  make build-run - Build and run"
	@echo "  make clean     - Remove build files"
