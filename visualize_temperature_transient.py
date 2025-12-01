#!/usr/bin/env python3
"""
Transient Temperature Distribution Visualization
Reads results from CSV and creates 2D visualization using matplotlib with time slider
Shows temperature evolution over time
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
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

def visualize_temperature_transient(csv_filename='results_transient.csv', grid_filename='grids/Test2_4_4_MixGrid.txt'):
    """
    Visualize transient temperature distribution from FEM results with time slider
    
    Args:
        csv_filename: Path to CSV file with node data
        grid_filename: Path to grid file with element connectivity
    """
    try:
        # Read the CSV file
        data = []
        with open(csv_filename, 'r') as f:
            lines = f.readlines()
            header = lines[0].strip().split(',')
            
            for line in lines[1:]:  # Skip header
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    node_id = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    temps = [float(t) for t in parts[3:]]  # All time steps
                    data.append([node_id, x, y] + temps)
        
        data = np.array(data)
        node_ids = data[:, 0].astype(int)
        x = data[:, 1]
        y = data[:, 2]
        temperatures = data[:, 3:]  # All time steps
        
        num_steps = temperatures.shape[1]
        
        print(f"\nLoaded {len(node_ids)} nodes from {csv_filename}")
        print(f"Number of time steps: {num_steps}")
        print(f"Temperature range: {temperatures.min():.2f}°C to {temperatures.max():.2f}°C\n")
        
        # Read element connectivity and boundary conditions
        elements, boundary_nodes = read_grid_file(grid_filename)
        print(f"Loaded {len(elements)} elements")
        print(f"Boundary nodes: {sorted(boundary_nodes)}\n")
        
        # Create figure with space for sliders
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # --- LAYOUT ADJUSTMENT ---
        # Reduced bottom margin to 0.25 (was 0.35) to reduce empty space
        plt.subplots_adjust(bottom=0.25, top=0.95)
        
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
        
        # Global temperature range for consistent colorbar
        temp_min = temperatures.min()
        temp_max = temperatures.max()
        
        # Storage for plot elements that need updating
        edge_lines = []
        
        # Draw element edges (lines between connected nodes)
        for element in elements:
            element_lines = []
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
                    line, = ax.plot(x_coords, y_coords, 'r-', linewidth=2, zorder=3, alpha=0.8)
                else:
                    line, = ax.plot(x_coords, y_coords, 'b-', linewidth=1.5, zorder=2, alpha=0.6)
                element_lines.append(line)
            edge_lines.append(element_lines)
        
        # Initial temperature (first time step)
        temp_initial = temperatures[:, 0]
        
        # Initial plot parameters
        initial_point_size = 500
        initial_badge_size = 9  # font size
        
        # Plot nodes with temperature-based colors
        scatter = ax.scatter(x, y, c=temp_initial, s=initial_point_size, cmap='hot', 
                           edgecolors='black', linewidth=2.5, zorder=5, 
                           vmin=temp_min, vmax=temp_max)
        
        # Add node labels with temperature values
        annotations = []
        for i, node_id in enumerate(node_ids):
            # Determine text color based on temperature
            temp_normalized = (temp_initial[i] - temp_min) / (temp_max - temp_min + 1e-10)
            
            if temp_normalized > 0.5:
                bg_color = 'darkred'
                text_color = 'white'
            else:
                bg_color = 'lightgray'
                text_color = 'black'
            
            # Position label above the node
            y_offset = (y_max - y_min) * 0.05
            
            ann = ax.annotate(f'Node {int(node_id)}\n{temp_initial[i]:.1f}°C', 
                       (x[i], y[i] + y_offset), 
                       ha='center', va='bottom',
                       fontsize=initial_badge_size, fontweight='bold',
                       color=text_color,
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor=bg_color,
                                edgecolor='black', 
                                alpha=1.0))
            annotations.append(ann)
        
        ax.set_xlabel('X coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_title('Temperature Distribution at Nodes - Time Step: 0 (t = 0.0 s)', 
                    fontsize=14, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8, color='gray')
        ax.set_aspect('equal')
        
        # Set plot limits with 20% padding
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, label='Temperature [°C]')
        cbar.ax.tick_params(labelsize=10)
        
        # --- SLIDER POSITIONS ---
        # Format: [left, bottom, width, height]
        # Width reduced to 0.5 (was 0.7)
        # Left set to 0.25 to center them ((1.0 - 0.5)/2)
        # Vertical positions adjusted for the tighter bottom margin
        
        # Point size slider (Top)
        ax_point_slider = plt.axes([0.25, 0.15, 0.5, 0.03], facecolor='lightgray')
        point_slider = Slider(ax_point_slider, 'Point Size', 50, 1000, valinit=initial_point_size, valstep=10)
        
        # Badge size slider (Middle)
        ax_badge_slider = plt.axes([0.25, 0.11, 0.5, 0.03], facecolor='lightgray')
        badge_slider = Slider(ax_badge_slider, 'Badge Size', 0, 20, valinit=initial_badge_size, valstep=1)
        
        # Time step slider (Bottom)
        ax_time_slider = plt.axes([0.25, 0.07, 0.5, 0.03], facecolor='lightgray')
        time_slider = Slider(ax_time_slider, 'Time Step', 0, num_steps - 1, valinit=0, valstep=1)
        
        # Update function for time slider
        def update_time(val):
            step = int(time_slider.val)
            temp_current = temperatures[:, step]
            
            # Update scatter plot colors
            scatter.set_array(temp_current)
            
            # Get current badge size
            badge_size = badge_slider.val
            
            # Update annotations
            for i, node_id in enumerate(node_ids):
                if badge_size == 0:
                    # Hide badges
                    annotations[i].set_visible(False)
                else:
                    # Show badges
                    annotations[i].set_visible(True)
                    
                    # Determine text color based on temperature
                    temp_normalized = (temp_current[i] - temp_min) / (temp_max - temp_min + 1e-10)
                    
                    if temp_normalized > 0.5:
                        bg_color = 'darkred'
                        text_color = 'white'
                    else:
                        bg_color = 'lightgray'
                        text_color = 'black'
                    
                    # Update annotation text and color
                    annotations[i].set_text(f'Node {int(node_id)}\n{temp_current[i]:.1f}°C')
                    annotations[i].set_color(text_color)
                    annotations[i].get_bbox_patch().set_facecolor(bg_color)
            
            # Update title with current time
            # Assume time step info from header or calculate
            time_value = step * 50.0  # Assuming 50s time step
            ax.set_title(f'Temperature Distribution at Nodes - Time Step: {step} (t = {time_value:.1f} s)', 
                        fontsize=14, fontweight='bold', pad=15)
            
            fig.canvas.draw_idle()
        
        # Update function for point size slider
        def update_point_size(val):
            point_size = point_slider.val
            scatter.set_sizes([point_size] * len(x))
            # Scale the border linewidth proportionally to point size
            # Base linewidth is 2.5 for size 500
            new_linewidth = 2.5 * (point_size / 500.0)
            scatter.set_linewidths(new_linewidth)
            fig.canvas.draw_idle()
        
        # Update function for badge size slider
        def update_badge_size(val):
            badge_size = badge_slider.val
            
            if badge_size == 0:
                # Hide all badges
                for ann in annotations:
                    ann.set_visible(False)
            else:
                # Show all badges and update font size
                for ann in annotations:
                    ann.set_visible(True)
                    ann.set_fontsize(badge_size)
            
            fig.canvas.draw_idle()
        
        time_slider.on_changed(update_time)
        point_slider.on_changed(update_point_size)
        badge_slider.on_changed(update_badge_size)
        
        # Show plot (opens GUI window)
        print("Opening visualization window...")
        print("Use the sliders to:")
        print("  - Time Step: Navigate through time steps")
        print("  - Point Size: Adjust the size of nodes (50-1000)")
        print("  - Badge Size: Adjust label size (0 = hide labels)")
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
        csv_file = 'results_transient.csv'
    
    # Try to find the grid file
    grid_file = 'grids/Test2_4_4_MixGrid.txt'
    if len(sys.argv) > 2:
        grid_file = sys.argv[2]
    
    visualize_temperature_transient(csv_file, grid_file)#!/usr/bin/env python3
"""
Transient Temperature Distribution Visualization
Reads results from CSV and creates 2D visualization using matplotlib with time slider
Shows temperature evolution over time
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
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

def visualize_temperature_transient(csv_filename='results_transient.csv', grid_filename='grids/Test2_4_4_MixGrid.txt'):
    """
    Visualize transient temperature distribution from FEM results with time slider
    
    Args:
        csv_filename: Path to CSV file with node data
        grid_filename: Path to grid file with element connectivity
    """
    try:
        # Read the CSV file
        data = []
        with open(csv_filename, 'r') as f:
            lines = f.readlines()
            header = lines[0].strip().split(',')
            
            for line in lines[1:]:  # Skip header
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    node_id = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    temps = [float(t) for t in parts[3:]]  # All time steps
                    data.append([node_id, x, y] + temps)
        
        data = np.array(data)
        node_ids = data[:, 0].astype(int)
        x = data[:, 1]
        y = data[:, 2]
        temperatures = data[:, 3:]  # All time steps
        
        num_steps = temperatures.shape[1]
        
        print(f"\nLoaded {len(node_ids)} nodes from {csv_filename}")
        print(f"Number of time steps: {num_steps}")
        print(f"Temperature range: {temperatures.min():.2f}°C to {temperatures.max():.2f}°C\n")
        
        # Read element connectivity and boundary conditions
        elements, boundary_nodes = read_grid_file(grid_filename)
        print(f"Loaded {len(elements)} elements")
        print(f"Boundary nodes: {sorted(boundary_nodes)}\n")
        
        # Create figure with space for sliders
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # --- LAYOUT ADJUSTMENT ---
        # Reduced bottom margin further to 0.20 to tighten the interface
        plt.subplots_adjust(bottom=0.20, top=0.95)
        
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
        
        # Global temperature range for consistent colorbar
        temp_min = temperatures.min()
        temp_max = temperatures.max()
        
        # Storage for plot elements that need updating
        edge_lines = []
        
        # Draw element edges (lines between connected nodes)
        for element in elements:
            element_lines = []
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
                    line, = ax.plot(x_coords, y_coords, 'r-', linewidth=2, zorder=3, alpha=0.8)
                else:
                    line, = ax.plot(x_coords, y_coords, 'b-', linewidth=1.5, zorder=2, alpha=0.6)
                element_lines.append(line)
            edge_lines.append(element_lines)
        
        # Initial temperature (first time step)
        temp_initial = temperatures[:, 0]
        
        # Initial plot parameters
        initial_point_size = 500
        initial_badge_size = 9  # font size
        
        # Plot nodes with temperature-based colors
        scatter = ax.scatter(x, y, c=temp_initial, s=initial_point_size, cmap='hot', 
                           edgecolors='black', linewidth=2.5, zorder=5, 
                           vmin=temp_min, vmax=temp_max)
        
        # Add node labels with temperature values
        annotations = []
        for i, node_id in enumerate(node_ids):
            # Determine text color based on temperature
            temp_normalized = (temp_initial[i] - temp_min) / (temp_max - temp_min + 1e-10)
            
            if temp_normalized > 0.5:
                bg_color = 'darkred'
                text_color = 'white'
            else:
                bg_color = 'lightgray'
                text_color = 'black'
            
            # Position label above the node
            y_offset = (y_max - y_min) * 0.05
            
            ann = ax.annotate(f'Node {int(node_id)}\n{temp_initial[i]:.1f}°C', 
                       (x[i], y[i] + y_offset), 
                       ha='center', va='bottom',
                       fontsize=initial_badge_size, fontweight='bold',
                       color=text_color,
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor=bg_color,
                                edgecolor='black', 
                                alpha=1.0))
            annotations.append(ann)
        
        ax.set_xlabel('X coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y coordinate [m]', fontsize=12, fontweight='bold')
        ax.set_title('Temperature Distribution at Nodes - Time Step: 0 (t = 0.0 s)', 
                    fontsize=14, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8, color='gray')
        ax.set_aspect('equal')
        
        # Set plot limits with 20% padding
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, label='Temperature [°C]')
        cbar.ax.tick_params(labelsize=10)
        
        # --- SLIDER POSITIONS ---
        # Format: [left, bottom, width, height]
        # Width: 0.5 (centered with left=0.25)
        # Height: 0.03
        # Vertical spacing reduced significantly (0.04 step between sliders)
        
        # Point size slider (Top: y=0.12)
        ax_point_slider = plt.axes([0.25, 0.12, 0.5, 0.03], facecolor='lightgray')
        point_slider = Slider(ax_point_slider, 'Point Size', 50, 1000, valinit=initial_point_size, valstep=10)
        
        # Badge size slider (Middle: y=0.08)
        ax_badge_slider = plt.axes([0.25, 0.08, 0.5, 0.03], facecolor='lightgray')
        badge_slider = Slider(ax_badge_slider, 'Badge Size', 0, 20, valinit=initial_badge_size, valstep=1)
        
        # Time step slider (Bottom: y=0.04)
        ax_time_slider = plt.axes([0.25, 0.04, 0.5, 0.03], facecolor='lightgray')
        time_slider = Slider(ax_time_slider, 'Time Step', 0, num_steps - 1, valinit=0, valstep=1)
        
        # Update function for time slider
        def update_time(val):
            step = int(time_slider.val)
            temp_current = temperatures[:, step]
            
            # Update scatter plot colors
            scatter.set_array(temp_current)
            
            # Get current badge size
            badge_size = badge_slider.val
            
            # Update annotations
            for i, node_id in enumerate(node_ids):
                if badge_size == 0:
                    # Hide badges
                    annotations[i].set_visible(False)
                else:
                    # Show badges
                    annotations[i].set_visible(True)
                    
                    # Determine text color based on temperature
                    temp_normalized = (temp_current[i] - temp_min) / (temp_max - temp_min + 1e-10)
                    
                    if temp_normalized > 0.5:
                        bg_color = 'darkred'
                        text_color = 'white'
                    else:
                        bg_color = 'lightgray'
                        text_color = 'black'
                    
                    # Update annotation text and color
                    annotations[i].set_text(f'Node {int(node_id)}\n{temp_current[i]:.1f}°C')
                    annotations[i].set_color(text_color)
                    annotations[i].get_bbox_patch().set_facecolor(bg_color)
            
            # Update title with current time
            # Assume time step info from header or calculate
            time_value = step * 50.0  # Assuming 50s time step
            ax.set_title(f'Temperature Distribution at Nodes - Time Step: {step} (t = {time_value:.1f} s)', 
                        fontsize=14, fontweight='bold', pad=15)
            
            fig.canvas.draw_idle()
        
        # Update function for point size slider
        def update_point_size(val):
            point_size = point_slider.val
            scatter.set_sizes([point_size] * len(x))
            # Scale the border linewidth proportionally to point size
            # Base linewidth is 2.5 for size 500
            new_linewidth = 2.5 * (point_size / 500.0)
            scatter.set_linewidths(new_linewidth)
            fig.canvas.draw_idle()
        
        # Update function for badge size slider
        def update_badge_size(val):
            badge_size = badge_slider.val
            
            if badge_size == 0:
                # Hide all badges
                for ann in annotations:
                    ann.set_visible(False)
            else:
                # Show all badges and update font size
                for ann in annotations:
                    ann.set_visible(True)
                    ann.set_fontsize(badge_size)
            
            fig.canvas.draw_idle()
        
        time_slider.on_changed(update_time)
        point_slider.on_changed(update_point_size)
        badge_slider.on_changed(update_badge_size)
        
        # Show plot (opens GUI window)
        print("Opening visualization window...")
        print("Use the sliders to:")
        print("  - Time Step: Navigate through time steps")
        print("  - Point Size: Adjust the size of nodes (50-1000)")
        print("  - Badge Size: Adjust label size (0 = hide labels)")
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
        csv_file = 'results_transient.csv'
    
    # Try to find the grid file
    grid_file = 'grids/Test2_4_4_MixGrid.txt'
    if len(sys.argv) > 2:
        grid_file = sys.argv[2]
    
    visualize_temperature_transient(csv_file, grid_file)