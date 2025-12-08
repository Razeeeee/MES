#!/usr/bin/env python3
"""
Generate plastic frames for window pane - add top and bottom frames
Keep existing window pane unchanged
"""

# Current window pane: 9cm wide (0.09m) x 80cm tall (0.80m)
# Window grid: 19 nodes wide x 81 nodes tall = 1539 nodes
# Window elements: 18 x 80 = 1440 elements

# Plastic frames: 15cm wide (0.15m), thickness 2cm (0.02m)
window_width = 0.09
window_height = 0.80
dx = 0.005  # 0.5cm horizontal spacing
dy = 0.01   # 1cm vertical spacing

frame_width = 0.15
frame_thickness = 0.02

# Frame positioning - centered horizontally
frame_start_x = (window_width - frame_width) / 2  # -3cm
frame_end_x = frame_start_x + frame_width  # 12cm

# Node counts
window_nx = 19  # 0 to 9cm
window_ny = 81  # 0 to 80cm
window_nodes = window_nx * window_ny  # 1539

frame_nx = int(frame_width / dx) + 1  # 31 nodes (0 to 15cm)
frame_ny = int(frame_thickness / dy) + 1  # 3 nodes (2cm)

bottom_frame_nodes = frame_nx * frame_ny  # 93
top_frame_nodes = frame_nx * frame_ny  # 93
total_nodes = window_nodes + bottom_frame_nodes + top_frame_nodes  # 1725

# Elements
window_elements = (window_nx - 1) * (window_ny - 1)  # 1440
frame_elements_per_section = (frame_nx - 1) * (frame_ny - 1)  # 60
total_elements = window_elements + 2 * frame_elements_per_section  # 1560

print(f"Total nodes: {total_nodes}")
print(f"Total elements: {total_elements}")

# Read existing file
with open('grids/WindowPane_Glass_Air_Glass.txt', 'r') as f:
    lines = f.readlines()

# Update header
for i, line in enumerate(lines):
    if line.startswith('Nodes number'):
        lines[i] = f'Nodes number {total_nodes}\n'
    elif line.startswith('Elements number'):
        lines[i] = f'Elements number {total_elements}\n'

# Find material section and add plastic material
material_idx = None
for i, line in enumerate(lines):
    if line.strip() == '*Material':
        material_idx = i
        break

# Insert plastic material after existing materials
insert_idx = material_idx + 3  # After Glass and Air
lines.insert(insert_idx, '3, Plastic, 0.2, 1200, 1500\n')

# Find where nodes end (before *Element line)
elem_section_idx = None
for i, line in enumerate(lines):
    if line.strip().startswith('*Element'):
        elem_section_idx = i
        break

# Generate bottom frame nodes (y from -0.02 to 0)
bottom_frame_node_list = []
node_id = window_nodes + 1
for j in range(frame_ny):
    y = -frame_thickness + j * dy
    for i in range(frame_nx):
        x = frame_start_x + i * dx
        bottom_frame_node_list.append(f"{node_id:7d}, {x:11.8f}, {y:11.8f}\n")
        node_id += 1

# Generate top frame nodes (y from 0.80 to 0.82)
top_frame_node_list = []
for j in range(frame_ny):
    y = window_height + j * dy
    for i in range(frame_nx):
        x = frame_start_x + i * dx
        top_frame_node_list.append(f"{node_id:7d}, {x:11.8f}, {y:11.8f}\n")
        node_id += 1

# Insert new nodes before *Element section
lines = lines[:elem_section_idx] + bottom_frame_node_list + top_frame_node_list + lines[elem_section_idx:]

# Find BC section
bc_idx = None
for i, line in enumerate(lines):
    if line.strip() == '*BC':
        bc_idx = i
        break

# Generate bottom frame elements
bottom_frame_elements_list = []
elem_id = window_elements + 1
base_node_id = window_nodes + 1

for j in range(frame_ny - 1):
    for i in range(frame_nx - 1):
        n1 = base_node_id + j * frame_nx + i
        n2 = n1 + 1
        n3 = n1 + frame_nx + 1
        n4 = n1 + frame_nx
        bottom_frame_elements_list.append(f"{elem_id:2d}, {n1:4d}, {n2:4d}, {n3:4d}, {n4:4d}, 3\n")
        elem_id += 1

# Generate top frame elements
top_frame_elements_list = []
base_node_id = window_nodes + bottom_frame_nodes + 1

for j in range(frame_ny - 1):
    for i in range(frame_nx - 1):
        n1 = base_node_id + j * frame_nx + i
        n2 = n1 + 1
        n3 = n1 + frame_nx + 1
        n4 = n1 + frame_nx
        top_frame_elements_list.append(f"{elem_id:2d}, {n1:4d}, {n2:4d}, {n3:4d}, {n4:4d}, 3\n")
        elem_id += 1

# Insert elements before *BC
lines = lines[:bc_idx] + bottom_frame_elements_list + top_frame_elements_list + lines[bc_idx:]

# Generate boundary conditions for frames
# Bottom frame: left edge, right edge, bottom edge, and top edge (touching window)
# Top frame: left edge, right edge, top edge, and bottom edge (touching window)
bc_nodes = []

# Bottom frame BC:
# - Bottom edge (j=0): all nodes
for i in range(frame_nx):
    node_id = window_nodes + i + 1
    bc_nodes.append(f"{node_id}:0")

# - Left edge (i=0): nodes except bottom corner (already added)
for j in range(1, frame_ny):
    node_id = window_nodes + j * frame_nx + 1
    bc_nodes.append(f"{node_id}:0")

# - Right edge (i=frame_nx-1): nodes except bottom corner
for j in range(1, frame_ny):
    node_id = window_nodes + j * frame_nx + frame_nx
    bc_nodes.append(f"{node_id}:20")

# - Top edge of bottom frame (j=frame_ny-1, touching window bottom)
for i in range(1, frame_nx - 1):
    node_id = window_nodes + (frame_ny - 1) * frame_nx + i + 1
    bc_nodes.append(f"{node_id}:10")

# Top frame BC:
base_top = window_nodes + bottom_frame_nodes

# - Bottom edge of top frame (j=0, touching window top)
for i in range(frame_nx):
    node_id = base_top + i + 1
    bc_nodes.append(f"{node_id}:10")

# - Left edge (i=0): nodes except already added
for j in range(1, frame_ny):
    node_id = base_top + j * frame_nx + 1
    bc_nodes.append(f"{node_id}:0")

# - Right edge (i=frame_nx-1): nodes except already added
for j in range(1, frame_ny):
    node_id = base_top + j * frame_nx + frame_nx
    bc_nodes.append(f"{node_id}:20")

# - Top edge (j=frame_ny-1): interior nodes only
for i in range(1, frame_nx - 1):
    node_id = base_top + (frame_ny - 1) * frame_nx + i + 1
    bc_nodes.append(f"{node_id}:0")

# Find and update BC section
bc_line_idx = None
for i, line in enumerate(lines):
    if line.strip() == '*BC':
        bc_line_idx = i
        break

# Get existing BC line
existing_bc = lines[bc_line_idx + 1].strip()

# Append new BC nodes
new_bc = existing_bc + ', ' + ', '.join(bc_nodes) + '\n'
lines[bc_line_idx + 1] = new_bc

# Write updated file
with open('grids/WindowPane_Glass_Air_Glass.txt', 'w') as f:
    f.writelines(lines)

print(f"✓ Added bottom frame: {bottom_frame_nodes} nodes, {frame_elements_per_section} elements")
print(f"✓ Added top frame: {top_frame_nodes} nodes, {frame_elements_per_section} elements")
print(f"✓ Total: {total_nodes} nodes, {total_elements} elements")
print(f"✓ Frame material (ID=3): Plastic (k=0.2, ρ=1200, c=1500)")
