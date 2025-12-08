#!/usr/bin/env python3
"""
Interactive Transient Temperature Visualization using Plotly
Advanced FEM heat transfer visualization with full interactivity:
- Pan and zoom
- Hover to see exact temperature and node info
- Toggle elements (nodes, labels, boundaries, mesh, contour plot)
- Time animation slider
- Export as standalone HTML

Compatible with Python 3.10+
Requires: plotly, numpy, pandas, scipy
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import griddata
import sys
import os


def load_transient_data(csv_file):
    """Load transient results from CSV file"""
    if not os.path.exists(csv_file):
        print(f"Error: File '{csv_file}' not found.")
        sys.exit(1)
    
    # Read using pandas for easier handling
    df = pd.read_csv(csv_file)
    
    # Extract basic info
    node_ids = df['node_id'].values
    x = df['x'].values
    y = df['y'].values
    
    # Extract time steps from column names
    time_cols = [col for col in df.columns if col.startswith('t_')]
    time_steps = np.array([float(col.split('_')[1]) for col in time_cols])
    
    # Extract temperature data
    temperatures = df[time_cols].values
    
    return node_ids, x, y, temperatures, time_steps


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
                element_nodes = [int(parts[i]) for i in range(1, 5)]
                elements.append(element_nodes)
        
        if in_bc_section and line:
            parts = [p.strip() for p in line.split(',')]
            for p in parts:
                if p:
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
            grid_files.sort(key=lambda f: os.path.getmtime(os.path.join(grid_dir, f)), reverse=True)
            return os.path.join(grid_dir, grid_files[0])
    return None


def create_mesh_edges(elements, node_ids, x, y, boundary_nodes):
    """Create edge lines for element visualization"""
    internal_edges_x = []
    internal_edges_y = []
    boundary_edges_x = []
    boundary_edges_y = []
    
    for element in elements:
        for i in range(4):
            node1_id = element[i]
            node2_id = element[(i + 1) % 4]
            
            idx1 = np.where(node_ids == node1_id)[0]
            idx2 = np.where(node_ids == node2_id)[0]
            
            if len(idx1) > 0 and len(idx2) > 0:
                idx1, idx2 = idx1[0], idx2[0]
                
                is_boundary = (node1_id in boundary_nodes) and (node2_id in boundary_nodes)
                
                if is_boundary:
                    boundary_edges_x.extend([x[idx1], x[idx2], None])
                    boundary_edges_y.extend([y[idx1], y[idx2], None])
                else:
                    internal_edges_x.extend([x[idx1], x[idx2], None])
                    internal_edges_y.extend([y[idx1], y[idx2], None])
    
    return internal_edges_x, internal_edges_y, boundary_edges_x, boundary_edges_y


def visualize_interactive(csv_file='data/transient_results.csv'):
    """
    Create interactive Plotly visualization with all features
    """
    print(f"\nLoading data from {csv_file}...")
    node_ids, x, y, temperatures, time_steps = load_transient_data(csv_file)
    
    # Load grid file for element connectivity
    grid_file = find_grid_file()
    elements, boundary_nodes = load_grid_file(grid_file) if grid_file else ([], set())
    
    num_nodes = len(node_ids)
    num_steps = len(time_steps)
    
    print(f"Loaded {num_nodes} nodes")
    print(f"Time steps: {num_steps}")
    print(f"Temperature range: {temperatures.min():.2f}°C to {temperatures.max():.2f}°C")
    print(f"Time range: {time_steps[0]:.1f}s to {time_steps[-1]:.1f}s")
    if grid_file:
        print(f"Grid file: {grid_file}")
        print(f"Elements: {len(elements)}")
        print(f"Boundary nodes: {len(boundary_nodes)}\n")
    
    # Global temperature range for consistent colorbar
    temp_min = temperatures.min()
    temp_max = temperatures.max()
    
    # Create figure
    fig = go.Figure()
    
    # Create mesh edges
    if elements:
        internal_x, internal_y, boundary_x, boundary_y = create_mesh_edges(
            elements, node_ids, x, y, boundary_nodes
        )
        
        # Add internal edges (thin, blue)
        if internal_x:
            fig.add_trace(go.Scatter(
                x=internal_x,
                y=internal_y,
                mode='lines',
                line=dict(color='rgba(100, 150, 200, 0.4)', width=1),
                hoverinfo='skip',
                showlegend=True,
                name='Internal Mesh',
                visible=True
            ))
        
        # Add boundary edges (thick, red)
        if boundary_x:
            fig.add_trace(go.Scatter(
                x=boundary_x,
                y=boundary_y,
                mode='lines',
                line=dict(color='rgba(255, 0, 0, 0.8)', width=3),
                hoverinfo='skip',
                showlegend=True,
                name='Boundary',
                visible=True
            ))
    
    # Determine marker sizes based on number of nodes
    if num_nodes < 20:
        marker_size = 30
        label_font_size = 10
    elif num_nodes < 100:
        marker_size = 20
        label_font_size = 9
    else:
        marker_size = 12
        label_font_size = 8
    
    # Create grid for contour interpolation
    # Create a finer grid for smooth contours
    grid_resolution = 100
    xi = np.linspace(x.min(), x.max(), grid_resolution)
    yi = np.linspace(y.min(), y.max(), grid_resolution)
    xi_mesh, yi_mesh = np.meshgrid(xi, yi)
    
    # Interpolate initial temperature for contour
    from scipy.interpolate import griddata
    initial_temp = temperatures[:, 0]
    zi_initial = griddata((x, y), initial_temp, (xi_mesh, yi_mesh), method='cubic')
    
    initial_time = time_steps[0]
    
    # Add contour plot (hidden by default)
    fig.add_trace(go.Contour(
        x=xi,
        y=yi,
        z=zi_initial,
        colorscale='RdYlBu_r',
        contours=dict(
            start=temp_min,
            end=temp_max,
            size=(temp_max - temp_min) / 20
        ),
        colorbar=dict(
            title="Temp [°C]",
            thickness=20,
            len=0.7,
            x=1.02,
            y=0.5
        ),
        hovertemplate='X: %{x:.4f}<br>Y: %{y:.4f}<br>Temp: %{z:.4f}°C<extra></extra>',
        showlegend=True,
        name='Contour Plot',
        visible='legendonly',  # Hidden by default
        showscale=False  # Don't show separate colorbar (will use nodes colorbar)
    ))
    
    # Create hover text for initial state
    hover_text = [
        f"Node {int(nid)}<br>" +
        f"Position: ({x[i]:.4f}, {y[i]:.4f})<br>" +
        f"Temperature: {temp:.4f}°C<br>" +
        f"Time: {initial_time:.1f}s" +
        (f"<br><b>BOUNDARY</b>" if int(nid) in boundary_nodes else "")
        for i, (nid, temp) in enumerate(zip(node_ids, initial_temp))
    ]
    
    # Add nodes with temperature colors (visible by default)
    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            size=marker_size,
            color=initial_temp,
            colorscale='RdYlBu_r',
            cmin=temp_min,
            cmax=temp_max,
            colorbar=dict(
                title="Temp [°C]",
                thickness=20,
                len=0.7,
                x=1.02
            ),
            line=dict(color='black', width=2)
        ),
        text=hover_text,
        hoverinfo='text',
        showlegend=True,
        name='Nodes',
        visible=True  # Visible by default, can be toggled
    ))
    
    # Add node labels
    label_text = [f"N{int(nid)}<br>{temp:.1f}°C" for nid, temp in zip(node_ids, initial_temp)]
    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='text',
        text=label_text,
        textposition='bottom center',
        textfont=dict(size=label_font_size, color='black', family='Arial Black'),
        hoverinfo='skip',
        showlegend=True,
        name='Labels',
        visible='legendonly'  # Hidden by default, can be toggled
    ))
    
    # Create frames for animation
    frames = []
    
    for step in range(num_steps):
        temp_current = temperatures[:, step]
        current_time = time_steps[step]
        
        # Interpolate temperature for contour at this time step
        zi_current = griddata((x, y), temp_current, (xi_mesh, yi_mesh), method='cubic')
        
        # Create hover text with detailed info
        hover_text = [
            f"Node {int(nid)}<br>" +
            f"Position: ({x[i]:.4f}, {y[i]:.4f})<br>" +
            f"Temperature: {temp:.4f}°C<br>" +
            f"Time: {current_time:.1f}s" +
            (f"<br><b>BOUNDARY</b>" if int(nid) in boundary_nodes else "")
            for i, (nid, temp) in enumerate(zip(node_ids, temp_current))
        ]
        
        # Node labels for this frame
        label_text = [f"N{int(nid)}<br>{temp:.1f}°C" for nid, temp in zip(node_ids, temp_current)]
        
        # Determine how many edge traces we have and where contour is
        num_edge_traces = len([t for t in fig.data if 'Mesh' in t.name or 'Boundary' in t.name])
        contour_trace_idx = num_edge_traces  # Contour is right after edges
        nodes_trace_idx = num_edge_traces + 1  # Nodes after contour
        labels_trace_idx = num_edge_traces + 2  # Labels after nodes
        
        frames.append(go.Frame(
            data=[
                go.Contour(z=zi_current),  # Update contour
                go.Scatter(
                    marker=dict(color=temp_current),
                    text=hover_text
                ),  # Update nodes
                go.Scatter(text=label_text)  # Update labels
            ],
            traces=[contour_trace_idx, nodes_trace_idx, labels_trace_idx],
            name=str(step),
            layout=go.Layout(
                title_text=f"Temperature Distribution - t = {current_time:.1f}s (Step {step + 1}/{num_steps})"
            )
        ))
    
    fig.frames = frames
    
    # Add animation controls
    fig.update_layout(
        updatemenus=[
            dict(
                type='buttons',
                showactive=False,
                x=0.05,
                y=0.02,
                xanchor='left',
                yanchor='bottom',
                buttons=[
                    dict(
                        label='▶ Play',
                        method='animate',
                        args=[None, dict(
                            frame=dict(duration=200, redraw=True),
                            fromcurrent=True,
                            mode='immediate',
                            transition=dict(duration=0)
                        )]
                    ),
                    dict(
                        label='⏸ Pause',
                        method='animate',
                        args=[[None], dict(
                            frame=dict(duration=0, redraw=False),
                            mode='immediate',
                            transition=dict(duration=0)
                        )]
                    )
                ]
            )
        ],
        sliders=[dict(
            active=0,
            yanchor='bottom',
            y=0.02,
            xanchor='left',
            x=0.25,
            currentvalue=dict(
                prefix='Time Step: ',
                visible=True,
                xanchor='left',
                font=dict(size=14)
            ),
            len=0.7,
            steps=[dict(
                args=[[f.name], dict(
                    frame=dict(duration=0, redraw=True),
                    mode='immediate',
                    transition=dict(duration=0)
                )],
                label=f"{time_steps[int(f.name)]:.0f}s",
                method='animate'
            ) for f in frames]
        )]
    )
    
    # Update layout
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    padding = 0.1
    
    fig.update_layout(
        title=dict(
            text=f"Temperature Distribution - t = {time_steps[0]:.1f}s (Step 1/{num_steps})",
            font=dict(size=18, family='Arial Black'),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=dict(text='X Position [m]', font=dict(size=14)),
            range=[x.min() - x_range * padding, x.max() + x_range * padding],
            scaleanchor='y',
            scaleratio=1,
            showgrid=True,
            gridcolor='lightgray',
            zeroline=False
        ),
        yaxis=dict(
            title=dict(text='Y Position [m]', font=dict(size=14)),
            range=[y.min() - y_range * padding, y.max() + y_range * padding],
            showgrid=True,
            gridcolor='lightgray',
            zeroline=False
        ),
        plot_bgcolor='#f5f5f5',
        paper_bgcolor='white',
        width=1200,
        height=900,
        hovermode='closest',
        legend=dict(
            orientation='v',
            yanchor='top',
            y=0.98,
            xanchor='left',
            x=0.02,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='black',
            borderwidth=1
        ),
        margin=dict(l=80, r=150, t=100, b=120)
    )
    
    print("Creating interactive visualization...")
    print("\nFeatures:")
    print("  • Pan: Click and drag")
    print("  • Zoom: Scroll or use zoom tools")
    print("  • Hover: See detailed node/temperature information")
    print("  • Toggle: Click legend items to show/hide layers:")
    print("    - Nodes (markers) - visible by default")
    print("    - Labels - hidden by default") 
    print("    - Contour Plot - hidden by default")
    print("    - Mesh/Boundary - visible by default")
    print("  • Animate: Use Play button or time slider")
    print("\nOpening browser...")
    
    # Save as HTML and open
    output_file = 'data/temperature_visualization.html'
    fig.write_html(output_file, auto_open=True)
    print(f"\n✓ Saved interactive visualization to: {output_file}")
    print("  You can share this HTML file - it's fully self-contained!\n")


if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'data/transient_results.csv'
    visualize_interactive(csv_file)
