#!/usr/bin/env python3
"""
Transient Temperature Visualization with Time Slider
Visualizes FEM heat transfer results over time with interactive controls
Compatible with Python 3.10+
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import sys
import os


def load_transient_data(csv_file):
    """Load transient results from CSV file"""
    if not os.path.exists(csv_file):
        print(f"Error: File '{csv_file}' not found.")
        sys.exit(1)
    
    # Read data
    data = []
    time_steps = []
    
    with open(csv_file, 'r') as f:
        header = f.readline().strip().split(',')
        
        # Extract time steps from header (skip node_id, x, y)
        for col in header[3:]:
            # Format: t_0.0, t_50.0, etc.
            time_val = float(col.split('_')[1])
            time_steps.append(time_val)
        
        # Read node data
        for line in f:
            parts = line.strip().split(',')
            node_id = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            temps = [float(t) for t in parts[3:]]
            data.append([node_id, x, y] + temps)
    
    data = np.array(data)
    node_ids = data[:, 0].astype(int)
    x = data[:, 1]
    y = data[:, 2]
    temperatures = data[:, 3:]  # All time steps
    
    return node_ids, x, y, temperatures, np.array(time_steps)


def load_grid_file(grid_file):
    """Load element connectivity and boundary conditions from grid file"""
    elements = []
    boundary_nodes = set()
    
    if not os.path.exists(grid_file):
        return elements, boundary_nodes
    
    with open(grid_file, 'r') as f:
        lines = f.readlines()
    
    in_element_section = False
    in_bc_section = False
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('*Element'):
            in_element_section = True
            in_bc_section = False
            continue
        elif line.startswith('*BC'):
            in_element_section = False
            in_bc_section = True
            continue
        elif line.startswith('*'):
            in_element_section = False
            in_bc_section = False
            continue
        
        if in_element_section and line:
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 5:
                # Element format: id, n1, n2, n3, n4, [materialId]
                element_nodes = [int(parts[i]) for i in range(1, 5)]
                elements.append(element_nodes)
        
        if in_bc_section and line:
            # Boundary condition nodes (supports both "nodeId" and "nodeId:temp" formats)
            parts = [p.strip() for p in line.split(',')]
            for p in parts:
                if p:
                    # Handle both old format (just nodeId) and new format (nodeId:temp)
                    if ':' in p:
                        node_id = int(p.split(':')[0])
                    else:
                        node_id = int(p)
                    boundary_nodes.add(node_id)
    
    return elements, boundary_nodes


def find_grid_file():
    """Try to find the grid file used for this simulation"""
    grid_dir = 'grids'
    if os.path.exists(grid_dir):
        grid_files = [f for f in os.listdir(grid_dir) if f.endswith('.txt')]
        if grid_files:
            # Return the most recently modified grid file
            grid_files.sort(key=lambda f: os.path.getmtime(os.path.join(grid_dir, f)), reverse=True)
            return os.path.join(grid_dir, grid_files[0])
    return None


def visualize_transient_temperature(csv_file='data/transient_results.csv'):
    """
    Create interactive visualization with time slider
    """
    print(f"\nLoading data from {csv_file}...")
    node_ids, x, y, temperatures, time_steps = load_transient_data(csv_file)
    
    # Load grid file for element connectivity
    grid_file = find_grid_file()
    elements, boundary_nodes = load_grid_file(grid_file) if grid_file else ([], set())
    
    num_nodes = len(node_ids)
    num_steps = temperatures.shape[1]
    
    print(f"Loaded {num_nodes} nodes")
    print(f"Time steps: {num_steps}")
    print(f"Temperature range: {temperatures.min():.2f}°C to {temperatures.max():.2f}°C")
    print(f"Time range: {time_steps[0]:.1f}s to {time_steps[-1]:.1f}s")
    if grid_file:
        print(f"Grid file: {grid_file}")
        print(f"Elements: {len(elements)}")
        print(f"Boundary nodes: {len(boundary_nodes)}\n")
    
    # Create figure with square aspect
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.subplots_adjust(bottom=0.12, top=0.95, left=0.08, right=0.92)
    
    # Set colors
    fig.patch.set_facecolor('white')
    ax.set_facecolor('#f5f5f5')
    
    # Calculate plot boundaries to make it square
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2
    
    # Use the larger range to make it square
    max_range = max(x_max - x_min, y_max - y_min)
    padding = max_range * 0.05
    square_half = (max_range + 2 * padding) / 2
    
    # Global temperature range for consistent colorbar
    temp_min = temperatures.min()
    temp_max = temperatures.max()
    
    # Draw element edges
    edge_lines = []
    has_boundary_edge = False
    if elements:
        for element in elements:
            # Draw edges: 0->1, 1->2, 2->3, 3->0
            for i in range(4):
                node1_id = element[i]
                node2_id = element[(i + 1) % 4]
                
                # Find node coordinates (node IDs are 1-based)
                idx1 = np.where(node_ids == node1_id)[0]
                idx2 = np.where(node_ids == node2_id)[0]
                
                if len(idx1) > 0 and len(idx2) > 0:
                    idx1, idx2 = idx1[0], idx2[0]
                    x_coords = [x[idx1], x[idx2]]
                    y_coords = [y[idx1], y[idx2]]
                    
                    # Check if this edge is on the boundary
                    is_boundary_edge = (node1_id in boundary_nodes) and (node2_id in boundary_nodes)
                    
                    # Draw line: red/thick for boundary, blue/thin for internal
                    if is_boundary_edge:
                        line, = ax.plot(x_coords, y_coords, 'r-', linewidth=3, 
                                       zorder=3, alpha=0.9, label='Boundary' if not has_boundary_edge else '')
                        has_boundary_edge = True
                    else:
                        line, = ax.plot(x_coords, y_coords, 'b-', linewidth=1.5, 
                                       zorder=2, alpha=0.5)
                    edge_lines.append(line)
    
    # Initial state (first time step)
    temp_initial = temperatures[:, 0]
    
    # Create scatter plot
    scatter = ax.scatter(x, y, c=temp_initial, s=500, cmap='RdYlBu_r',
                        edgecolors='black', linewidth=2.5, zorder=5,
                        vmin=temp_min, vmax=temp_max)
    
    # Add node labels
    annotations = []
    for i, node_id in enumerate(node_ids):
        ann = ax.annotate(
            f'{int(node_id)}\n{temp_initial[i]:.1f}°C',
            (x[i], y[i]),
            ha='center', va='top',
            fontsize=7, fontweight='bold',
            color='black',
            xytext=(0, -18),
            textcoords='offset points',
            bbox=dict(boxstyle='round,pad=0.3',
                     facecolor='white',
                     edgecolor='black',
                     alpha=1.0)
        )
        annotations.append(ann)
    
    # Axis settings
    ax.set_xlabel('X Position [m]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y Position [m]', fontsize=12, fontweight='bold')
    ax.set_title(f'Temperature Distribution - t = {time_steps[0]:.1f}s',
                fontsize=14, fontweight='bold', pad=15)
    
    # Set limits with minimal padding (horizontal stretch)
    x_padding = (x_max - x_min) * 0.05
    y_padding = (y_max - y_min) * 0.15
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding, y_max + y_padding)
    ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, color='gray')
    
    # Add legend if we have boundary edges
    if has_boundary_edge:
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='Temperature [°C]', shrink=0.8)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_ticks(np.linspace(temp_min, temp_max, 5))
    cbar.set_ticklabels([f'{t:.1f}' for t in np.linspace(temp_min, temp_max, 5)])
    
    # Time slider
    ax_slider = plt.axes([0.25, 0.05, 0.5, 0.03], facecolor='lightgray')
    time_slider = Slider(
        ax_slider,
        'Time',
        0,
        num_steps - 1,
        valinit=0,
        valstep=1,
        valfmt='%0.0f'
    )
    
    # Update function
    def update_time(val):
        step = int(time_slider.val)
        temp_current = temperatures[:, step]
        current_time = time_steps[step]
        
        # Update scatter colors
        scatter.set_array(temp_current)
        
        # Update annotations
        for i in range(len(annotations)):
            annotations[i].set_text(f'{int(node_ids[i])}\n{temp_current[i]:.1f}°C')
        
        # Update title
        ax.set_title(
            f'Temperature Distribution - t = {current_time:.1f}s (Step {step + 1}/{num_steps})',
            fontsize=14, fontweight='bold', pad=15
        )
        
        fig.canvas.draw_idle()
    
    time_slider.on_changed(update_time)
    
    print("Opening visualization window...")
    print("Use the slider to navigate through time steps.")
    print("Close the window to exit.\n")
    
    plt.show()


if __name__ == "__main__":
    # Get filename from command line or use default
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'data/transient_results.csv'
    visualize_transient_temperature(csv_file)
