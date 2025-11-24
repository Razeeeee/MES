#!/usr/bin/env python3
"""
Temperature Distribution Visualization
Reads results from CSV and creates 2D visualization using matplotlib
Shows only actual node positions with their temperatures
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
import os

def read_grid_file(grid_filename):
    """Read element connectivity and boundary conditions from grid file"""
    elements = []
    boundary_nodes = set()
    
    if not os.path.exists(grid_filename):
        return elements, boundary_nodes
    
    with open(grid_filename, 'r') as f:
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
                # Element format: id, node1, node2, node3, node4
                element_nodes = [int(parts[i]) for i in range(1, 5)]
                elements.append(element_nodes)
        
        if in_bc_section and line:
            # Boundary condition nodes
            parts = [p.strip() for p in line.split(',')]
            for p in parts:
                if p:
                    boundary_nodes.add(int(p))
    
    return elements, boundary_nodes

def visualize_temperature(csv_filename='results.csv', grid_filename='grids/Test2_4_4_MixGrid.txt'):
    """
    Visualize temperature distribution from FEM results - nodes only
    
    Args:
        csv_filename: Path to CSV file with node data
    """
    try:
        # Read the CSV file manually (no pandas needed)
        data = []
        with open(csv_filename, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip header
                parts = line.strip().split(',')
                if len(parts) == 4:
                    node_id = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    temp = float(parts[3])
                    data.append([node_id, x, y, temp])
        
        data = np.array(data)
        node_ids = data[:, 0].astype(int)
        x = data[:, 1]
        y = data[:, 2]
        temp = data[:, 3]
        
        print(f"\nLoaded {len(node_ids)} nodes from {csv_filename}")
        print(f"Temperature range: {temp.min():.2f}°C to {temp.max():.2f}°C\n")
        
        # Read element connectivity and boundary conditions
        elements, boundary_nodes = read_grid_file(grid_filename)
        print(f"Loaded {len(elements)} elements")
        print(f"Boundary nodes: {sorted(boundary_nodes)}\n")
        
        # Create single plot showing nodes with temperatures
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Set background color
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        
        # Calculate plot boundaries with 20% expansion
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        x_padding = x_range * 0.20
        y_padding = y_range * 0.20
        
        # Draw element edges (lines between connected nodes)
        for element in elements:
            # Each element has 4 nodes forming a quadrilateral
            # Draw edges: 0->1, 1->2, 2->3, 3->0
            for i in range(4):
                node1_id = element[i]
                node2_id = element[(i + 1) % 4]
                
                # Find node coordinates (node IDs are 1-based)
                idx1 = np.where(node_ids == node1_id)[0][0]
                idx2 = np.where(node_ids == node2_id)[0][0]
                
                x_coords = [x[idx1], x[idx2]]
                y_coords = [y[idx1], y[idx2]]
                
                # Determine if this edge is on the boundary
                is_boundary_edge = (node1_id in boundary_nodes) and (node2_id in boundary_nodes)
                
                # Draw line: red for boundary edges, blue for internal edges
                if is_boundary_edge:
                    ax.plot(x_coords, y_coords, 'r-', linewidth=2, zorder=3, alpha=0.8)
                else:
                    ax.plot(x_coords, y_coords, 'b-', linewidth=1.5, zorder=2, alpha=0.6)
        
        # Plot nodes with temperature-based colors
        scatter = ax.scatter(x, y, c=temp, s=500, cmap='hot', 
                           edgecolors='black', linewidth=2.5, zorder=5, 
                           vmin=temp.min(), vmax=temp.max())
        
        # Add node labels with temperature values (positioned above nodes)
        for i, node_id in enumerate(node_ids):
            # Determine text color based on temperature
            temp_normalized = (temp[i] - temp.min()) / (temp.max() - temp.min() + 1e-10)
            
            # Use contrasting colors: dark background with white text, or light background with black text
            if temp_normalized > 0.5:
                bg_color = 'darkred'
                text_color = 'white'
            else:
                bg_color = 'lightgray'
                text_color = 'black'
            
            # Position label above the node (increased offset)
            y_offset = (y_max - y_min) * 0.05  # Offset proportional to plot size
            
            ax.annotate(f'Node {int(node_id)}\n{temp[i]:.1f}°C', 
                       (x[i], y[i] + y_offset), 
                       ha='center', va='bottom',
                       fontsize=9, fontweight='bold',
                       color=text_color,
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor=bg_color,
                                edgecolor='black', 
                                alpha=1.0))
        
        ax.set_xlabel('X coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_title('Temperature Distribution at Nodes', fontsize=14, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8, color='gray')
        ax.set_aspect('equal')
        
        # Set plot limits with 20% padding
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, label='Temperature [°C]')
        cbar.ax.tick_params(labelsize=10)
        
        # Adjust layout
        plt.tight_layout()
        
        # Show plot (opens GUI window)
        print("Opening visualization window...")
        print("Close the window to continue...\n")
        plt.show()
        
    except FileNotFoundError:
        print(f"Error: File '{csv_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Get filename from command line argument or use default
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'results.csv'
    
    # Try to find the grid file
    grid_file = 'grids/Test2_4_4_MixGrid.txt'
    if len(sys.argv) > 2:
        grid_file = sys.argv[2]
    
    visualize_temperature(csv_file, grid_file)
