.PHONY: build run dev clean init docs visualize help

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
dev:
	@chmod +x scripts/build-run.sh
	@./scripts/build-run.sh

# Visualize transient results (legacy matplotlib version)
visualize:
	@if [ ! -f data/transient_results.csv ]; then \
		echo "Error: data/transient_results.csv not found. Run 'make dev' first."; \
		exit 1; \
	fi
	@echo "Checking Python dependencies..."
	@py -3.10 -c "import numpy, matplotlib" 2>/dev/null || \
		(echo "Error: Missing required Python packages." && \
		 echo "Install with: py -3.10 -m pip install numpy matplotlib" && \
		 exit 1)
	@echo "Launching transient temperature visualization..."
	@py -3.10 scripts/visualize_transient.py

# Interactive visualization with Plotly (recommended)
viz:
	@if [ ! -f data/transient_results.csv ]; then \
		echo "Error: data/transient_results.csv not found. Run 'make dev' first."; \
		exit 1; \
	fi
	@echo "Checking Python dependencies..."
	@py -3.10 -c "import numpy, pandas, plotly, scipy" 2>/dev/null || \
		(echo "Error: Missing required Python packages." && \
		 echo "Install with: py -3.10 -m pip install numpy pandas plotly scipy" && \
		 exit 1)
	@echo "Launching interactive visualization..."
	@py -3.10 scripts/visualize_interactive.py

# Generate window frame grid
generate-frame:
	@echo "Generating window system grid..."
	@py -3.10 scripts/generate_window_system.py

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
	@echo "  make init           - Set script permissions"
	@echo "  make build          - Build project"
	@echo "  make run            - Run executable"
	@echo "  make dev            - Build and run"
	@echo "  make generate-frame - Generate window frame grid"
	@echo "  make viz            - Interactive visualization (Plotly - RECOMMENDED)"
	@echo "  make visualize      - Legacy visualization (Matplotlib)"
	@echo "  make clean          - Remove build files"
	@echo "  make docs           - Compile LaTeX documentation"
