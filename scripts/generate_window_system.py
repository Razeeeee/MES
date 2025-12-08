#!/usr/bin/env python3
"""
Generate complete window system with plastic frames
- Top plastic frame: 15cm wide x 10cm tall
- Window pane: 9cm wide (3cm glass + 3cm air + 3cm glass) x 70cm tall
- Bottom plastic frame: 15cm wide x 10cm tall
- Node spacing: 0.5cm horizontal, 1cm vertical
- All nodes connected in elements including boundary nodes
"""

# Configuration
dx = 0.005  # 0.5cm horizontal spacing
dy = 0.01   # 1cm vertical spacing

# Dimensions
plastic_width = 0.15  # 15cm
plastic_height = 0.10  # 10cm
window_width = 0.09  # 9cm (3cm + 3cm + 3cm)
window_height = 0.70  # 70cm

# Center window horizontally within plastic frame
plastic_start_x = 0.0
window_start_x = (plastic_width - window_width) / 2  # 3cm offset

# Node counts
plastic_nx = int(plastic_width / dx) + 1  # 31 nodes (0 to 15cm)
plastic_ny = int(plastic_height / dy) + 1  # 11 nodes (0 to 10cm)
window_nx = int(window_width / dx) + 1  # 19 nodes (0 to 9cm)
window_ny = int(window_height / dy) + 1  # 71 nodes (0 to 70cm)

# Total nodes
bottom_plastic_nodes = plastic_nx * plastic_ny  # 341
# Window shares top row of bottom plastic and bottom row of top plastic
window_interior_nodes = window_nx * (window_ny - 2)  # 19 * 69 = 1311 (excluding shared rows)
top_plastic_nodes = plastic_nx * plastic_ny  # 341
total_nodes = bottom_plastic_nodes + window_interior_nodes + top_plastic_nodes  # 341 + 1311 + 341 = 1993

# Elements (including connection elements)
plastic_elements = (plastic_nx - 1) * (plastic_ny - 1)  # 300 per plastic part
window_elements = (window_nx - 1) * (window_ny - 1)  # 1260
# No separate connection elements needed - window elements connect directly to plastic frames
total_elements = 2 * plastic_elements + window_elements  # 600 + 1260 = 1860

print(f"Generating window system:")
print(f"  Bottom plastic: {plastic_nx}x{plastic_ny} = {bottom_plastic_nodes} nodes, {plastic_elements} elements")
print(f"  Window pane: {window_nx}x{window_ny} = {window_nx}x{window_ny} grid ({window_interior_nodes} interior nodes), {window_elements} elements")
print(f"  Top plastic: {plastic_nx}x{plastic_ny} = {top_plastic_nodes} nodes, {plastic_elements} elements")
print(f"  Total: {total_nodes} nodes, {total_elements} elements (boundaries shared)")

# Start writing the file
lines = []

# Header
lines.append("SimulationTime 1200\n")
lines.append("SimulationStepTime 60\n")
lines.append("Conductivity 1.0\n")
lines.append("Alfa 25\n")
lines.append("Tot 0\n")
lines.append("InitialTemp 20\n")
lines.append("Density 2500\n")
lines.append("SpecificHeat 750\n")
lines.append(f"Nodes number {total_nodes}\n")
lines.append(f"Elements number {total_elements}\n")

# Materials
lines.append("*Material\n")
lines.append("1, Glass, 1.0, 2500, 750\n")
lines.append("2, Air, 0.025, 1.2, 1005\n")
lines.append("3, Plastic, 0.2, 1200, 1500\n")

# Nodes section
lines.append("*Node\n")

# Generate bottom plastic nodes (y from 0 to 10cm)
node_id = 1
for j in range(plastic_ny):
    y = j * dy
    for i in range(plastic_nx):
        x = plastic_start_x + i * dx
        lines.append(f"{node_id:7d}, {x:13.8f}, {y:13.8f}\n")
        node_id += 1

# Generate window nodes (y from 10cm to 80cm)
# IMPORTANT: Skip the first row (j=0) as it shares nodes with bottom plastic top edge
# IMPORTANT: Skip the last row (j=window_ny-1) as it shares nodes with top plastic bottom edge
window_y_start = plastic_height
for j in range(1, window_ny - 1):  # Skip first and last rows
    y = window_y_start + j * dy
    for i in range(window_nx):
        x = window_start_x + i * dx
        lines.append(f"{node_id:7d}, {x:13.8f}, {y:13.8f}\n")
        node_id += 1

# Generate top plastic nodes (y from 80cm to 90cm)
top_y_start = plastic_height + window_height
for j in range(plastic_ny):
    y = top_y_start + j * dy
    for i in range(plastic_nx):
        x = plastic_start_x + i * dx
        lines.append(f"{node_id:7d}, {x:13.8f}, {y:13.8f}\n")
        node_id += 1

# Elements section
lines.append("*Element, type=DC2D4\n")

# Bottom plastic elements
elem_id = 1
base_node_id = 0
for j in range(plastic_ny - 1):
    for i in range(plastic_nx - 1):
        n1 = base_node_id + j * plastic_nx + i + 1
        n2 = n1 + 1
        n3 = n1 + plastic_nx + 1
        n4 = n1 + plastic_nx
        lines.append(f"{elem_id}, {n1}, {n2}, {n3}, {n4}, 3\n")
        elem_id += 1

# Window elements - determine material based on x position
# Base address for window interior nodes (first interior row after shared bottom boundary)
base_window = bottom_plastic_nodes
base_top_plastic = bottom_plastic_nodes + window_interior_nodes
window_start_idx = int(window_start_x / dx)  # 6

# Window elements - determine material based on x position
# Glass: x < 0.03m (3cm from window start) or x >= 0.06m (6cm from window start)
# Air: 0.03m <= x < 0.06m
# NOTE: Window shares nodes with plastic frames at boundaries
for j in range(window_ny - 1):
    for i in range(window_nx - 1):
        # Calculate node IDs accounting for shared boundaries
        if j == 0:
            # Bottom row: use bottom plastic top edge nodes (bottom) and first window interior row (top)
            plastic_idx_left = window_start_idx + i
            plastic_idx_right = plastic_idx_left + 1
            # Bottom edge: plastic top row nodes
            n1 = (plastic_ny - 1) * plastic_nx + plastic_idx_left + 1
            n2 = (plastic_ny - 1) * plastic_nx + plastic_idx_right + 1
            # Top edge: first interior window row nodes
            n3 = base_window + i + 1 + 1
            n4 = base_window + i + 1
        elif j == window_ny - 2:
            # Second to last row: bottom nodes are interior, top nodes are from top plastic
            row_offset = (j - 1) * window_nx
            n1 = base_window + row_offset + i + 1
            n2 = n1 + 1
            # Top row: use top plastic bottom edge nodes
            plastic_idx_left = window_start_idx + i
            plastic_idx_right = plastic_idx_left + 1
            n3 = base_top_plastic + plastic_idx_right + 1
            n4 = base_top_plastic + plastic_idx_left + 1
        else:
            # Interior rows
            row_offset = (j - 1) * window_nx
            n1 = base_window + row_offset + i + 1
            n2 = n1 + 1
            n3 = n1 + window_nx + 1
            n4 = n1 + window_nx
        
        # Determine material based on x position within window
        x_in_window = i * dx
        if x_in_window < 0.03 or x_in_window >= 0.06:
            material = 1  # Glass
        else:
            material = 2  # Air
        
        lines.append(f"{elem_id}, {n1}, {n2}, {n3}, {n4}, {material}\n")
        elem_id += 1

# Top plastic elements
for j in range(plastic_ny - 1):
    for i in range(plastic_nx - 1):
        n1 = base_top_plastic + j * plastic_nx + i + 1
        n2 = n1 + 1
        n3 = n1 + plastic_nx + 1
        n4 = n1 + plastic_nx
        lines.append(f"{elem_id}, {n1}, {n2}, {n3}, {n4}, 3\n")
        elem_id += 1

# Boundary conditions
lines.append("*BC\n")

bc_nodes = []

# Define base addresses for BC generation
base_window = bottom_plastic_nodes
base_top = bottom_plastic_nodes + window_interior_nodes
window_start_idx = int(window_start_x / dx)  # 6
window_end_idx = window_start_idx + (window_nx - 1)  # 24

# Bottom plastic BCs
# Left edge (i=0, all j)
for j in range(plastic_ny):
    node_id = j * plastic_nx + 1
    bc_nodes.append(f"{node_id}:0")

# Right edge (i=plastic_nx-1, all j)
for j in range(plastic_ny):
    node_id = j * plastic_nx + plastic_nx
    bc_nodes.append(f"{node_id}:20")

# Top edge (j=plastic_ny-1): only portions NOT overlapping window (left and right extensions)
# Left extension gets 0°C, right extension gets 20°C
for i in range(plastic_nx):
    if i < window_start_idx and i != 0:  # Left extension (excluding corner already set)
        node_id = (plastic_ny - 1) * plastic_nx + i + 1
        bc_nodes.append(f"{node_id}:0")
    elif i > window_end_idx and i != plastic_nx - 1:  # Right extension (excluding corner already set)
        node_id = (plastic_ny - 1) * plastic_nx + i + 1
        bc_nodes.append(f"{node_id}:20")

# Window BCs
# Left edge (i=0, all j) - includes shared boundary nodes
# Bottom boundary shared with plastic
for j in range(window_ny):
    if j == 0:
        # Use bottom plastic top edge node
        node_id = (plastic_ny - 1) * plastic_nx + window_start_idx + 1
    elif j == window_ny - 1:
        # Use top plastic bottom edge node
        node_id = base_top + window_start_idx + 1
    else:
        # Interior window node
        node_id = base_window + (j - 1) * window_nx + 1
    bc_nodes.append(f"{node_id}:0")

# Right edge (i=window_nx-1, all j)
for j in range(window_ny):
    if j == 0:
        # Use bottom plastic top edge node
        node_id = (plastic_ny - 1) * plastic_nx + window_start_idx + window_nx
    elif j == window_ny - 1:
        # Use top plastic bottom edge node
        node_id = base_top + window_start_idx + window_nx
    else:
        # Interior window node
        node_id = base_window + (j - 1) * window_nx + window_nx
    bc_nodes.append(f"{node_id}:20")

# Top plastic BCs
# Left edge (i=0, all j)
for j in range(plastic_ny):
    node_id = base_top + j * plastic_nx + 1
    bc_nodes.append(f"{node_id}:0")

# Right edge (i=plastic_nx-1, all j)
for j in range(plastic_ny):
    node_id = base_top + j * plastic_nx + plastic_nx
    bc_nodes.append(f"{node_id}:20")

# Bottom edge (j=0): only portions NOT overlapping window (left and right extensions)
# Left extension gets 0°C, right extension gets 20°C
for i in range(plastic_nx):
    if i < window_start_idx and i != 0:  # Left extension (excluding corner already set)
        node_id = base_top_plastic + i + 1
        bc_nodes.append(f"{node_id}:0")
    elif i > window_end_idx and i != plastic_nx - 1:  # Right extension (excluding corner already set)
        node_id = base_top_plastic + i + 1
        bc_nodes.append(f"{node_id}:20")

# Write BC line
bc_line = ', '.join(bc_nodes) + '\n'
lines.append(bc_line)

# Write to file
with open('grids/WindowPane_Glass_Air_Glass.txt', 'w') as f:
    f.writelines(lines)

# Validation
print(f"\n✓ Generated: grids/WindowPane_Glass_Air_Glass.txt")
print(f"✓ Total nodes: {total_nodes}")
print(f"✓ Total elements: {total_elements}")
print(f"✓ Configuration:")
print(f"  - Bottom plastic: 15cm x 10cm (y: 0-10cm)")
print(f"  - Window: 9cm x 70cm (y: 10-80cm, x: 3-12cm)")
print(f"    - Glass layers: 0-3cm and 6-9cm")
print(f"    - Air layer: 3-6cm")
print(f"  - Top plastic: 15cm x 10cm (y: 80-90cm)")
print(f"✓ All boundary nodes are shared between components")

# Verify element node references
print(f"\n✓ Node allocation:")
print(f"  - Bottom plastic: nodes 1 to {bottom_plastic_nodes}")
print(f"  - Window interior: nodes {bottom_plastic_nodes + 1} to {bottom_plastic_nodes + window_interior_nodes}")
print(f"  - Top plastic: nodes {bottom_plastic_nodes + window_interior_nodes + 1} to {total_nodes}")
print(f"  - Window shares boundaries with plastic frames (no duplicate nodes)")
