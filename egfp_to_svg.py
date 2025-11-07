#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pymol-open-source",
#     "matplotlib",
#     "numpy",
#     "biopython",
# ]
# ///

"""
Generate a 2D SVG representation of EGFP (PDB: 4EUL) beta barrel structure.

Usage:
    uv run egfp_to_svg.py

This will generate 'egfp_structure.svg' in the current directory.

This script loads the EGFP molecular structure from PDB, uses PyMOL for
visualization and orientation, then projects the 3D coordinates to 2D
and exports as an SVG vector graphic.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
import pymol
from pymol import cmd
from Bio.PDB import PDBParser, PDBIO


def fetch_structure(pdb_id: str = "4EUL"):
    """Fetch structure using PyMOL."""
    print(f"Fetching {pdb_id} from PDB...")
    cmd.fetch(pdb_id, name="egfp")
    cmd.remove("solvent")  # Remove water molecules
    print("Structure loaded successfully")


def setup_view():
    """Set up a nice view of the beta barrel."""
    print("Setting up view...")
    cmd.hide("everything")
    cmd.show("cartoon", "egfp")
    cmd.color("green", "egfp")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)

    # Orient to show the barrel nicely
    cmd.orient("egfp")
    cmd.turn("y", 90)  # Rotate to see the barrel from the side
    cmd.zoom("egfp", 2)


def get_coordinates_2d():
    """
    Extract 3D coordinates from PyMOL's current view and project to 2D.
    Returns CA atoms and secondary structure info.
    """
    print("Extracting coordinates...")

    # Get the model in current orientation
    model = cmd.get_model("egfp and name CA")  # C-alpha atoms only

    coords = []
    residues = []

    for atom in model.atom:
        coords.append([atom.coord[0], atom.coord[1], atom.coord[2]])
        residues.append({
            'index': atom.resi_number,
            'name': atom.resn,
            'ss': atom.ss,  # Secondary structure
        })

    coords = np.array(coords)

    # Get the current view matrix from PyMOL
    view = cmd.get_view()
    rotation_matrix = np.array(view[:9]).reshape(3, 3)

    # Apply rotation and project to 2D (just take x, y)
    rotated = coords @ rotation_matrix.T
    coords_2d = rotated[:, :2]  # Take x, y coordinates

    return coords_2d, residues


def create_svg(coords_2d, residues, output_file: str = "egfp_structure.svg"):
    """Create an SVG representation of the structure."""
    print(f"Creating SVG: {output_file}")

    fig, ax = plt.subplots(figsize=(12, 12), dpi=100)
    ax.set_aspect('equal')

    # Color map for secondary structure
    ss_colors = {
        'H': '#FF6B6B',  # Helix - red
        'S': '#4ECDC4',  # Sheet (beta strand) - cyan
        'L': '#95E1D3',  # Loop - light green
        '': '#AAAAAA',   # Unknown - gray
    }

    # Draw backbone as connected line segments with colors based on secondary structure
    segments = []
    colors = []

    for i in range(len(coords_2d) - 1):
        segment = [coords_2d[i], coords_2d[i + 1]]
        segments.append(segment)

        # Color based on secondary structure
        ss = residues[i].get('ss', '')
        colors.append(ss_colors.get(ss, '#AAAAAA'))

    # Create line collection for the backbone
    lc = LineCollection(segments, colors=colors, linewidths=3, alpha=0.8)
    ax.add_collection(lc)

    # Add CA atoms as small circles
    for i, (coord, res) in enumerate(zip(coords_2d, residues)):
        ss = res.get('ss', '')
        color = ss_colors.get(ss, '#AAAAAA')
        circle = Circle(coord, radius=0.3, color=color, alpha=0.6, zorder=2)
        ax.add_patch(circle)

    # Highlight beta strands (the barrel structure) with thicker lines
    beta_coords = [coords_2d[i] for i, res in enumerate(residues) if res.get('ss') == 'S']
    if beta_coords:
        beta_coords = np.array(beta_coords)
        ax.scatter(beta_coords[:, 0], beta_coords[:, 1],
                  c='#4ECDC4', s=50, alpha=0.9, zorder=3,
                  edgecolors='#2C3E50', linewidths=1,
                  label='β-strands (barrel)')

    # Add title and legend
    ax.set_title('EGFP Structure (PDB: 4EUL) - 2D Projection',
                fontsize=16, fontweight='bold', pad=20)

    legend_elements = [
        plt.Line2D([0], [0], color='#FF6B6B', lw=4, label='α-helix'),
        plt.Line2D([0], [0], color='#4ECDC4', lw=4, label='β-strand (barrel)'),
        plt.Line2D([0], [0], color='#95E1D3', lw=4, label='Loop/coil'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12)

    # Remove axes for cleaner look
    ax.set_xlim(coords_2d[:, 0].min() - 5, coords_2d[:, 0].max() + 5)
    ax.set_ylim(coords_2d[:, 1].min() - 5, coords_2d[:, 1].max() + 5)
    ax.axis('off')

    # Save as SVG
    plt.tight_layout()
    plt.savefig(output_file, format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print(f"✓ SVG saved to: {output_file}")

    # Also show stats
    n_helix = sum(1 for r in residues if r.get('ss') == 'H')
    n_sheet = sum(1 for r in residues if r.get('ss') == 'S')
    n_loop = len(residues) - n_helix - n_sheet

    print(f"\nStructure statistics:")
    print(f"  Total residues: {len(residues)}")
    print(f"  α-helices: {n_helix} residues")
    print(f"  β-strands: {n_sheet} residues (forming the barrel)")
    print(f"  Loops/coils: {n_loop} residues")

    plt.close()


def main() -> None:
    """Main function to generate EGFP SVG."""
    print("=" * 60)
    print("EGFP 3D Structure → 2D SVG Converter")
    print("=" * 60)
    print()

    # Initialize PyMOL in quiet mode (no GUI)
    pymol.finish_launching(['pymol', '-cq'])

    try:
        # Fetch and prepare structure
        fetch_structure("4EUL")
        setup_view()

        # Extract 2D coordinates
        coords_2d, residues = get_coordinates_2d()

        # Create SVG
        output_file = "egfp_structure.svg"
        create_svg(coords_2d, residues, output_file)

        print(f"\n{'=' * 60}")
        print("✓ Complete! Check out {output_file}")
        print(f"{'=' * 60}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        cmd.quit()


if __name__ == "__main__":
    main()
