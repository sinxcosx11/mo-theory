# mo_diagram_generator/plotter.py
import matplotlib.pyplot as plt
import numpy as np
from .orbitals import AtomicOrbital, MolecularOrbital
from typing import List, Dict, Tuple, Set
import math
import sys

# --- Plotting Configuration ---
AO_COLOR = 'cornflowerblue'
MO_COLOR = 'black'
# MODIFIED Colors/Width
BONDING_COLOR = 'darkblue' # Was 'darkgreen'
ANTIBONDING_COLOR = 'crimson' # Was 'firebrick'
NONBONDING_COLOR = 'dimgrey'
LINE_COLOR = 'darkgrey'
ELECTRON_COLOR = 'black'
AO_LINE_WIDTH = 1.5
MO_LINE_WIDTH = 2.0 # MODIFIED - Was 1.5
CONNECTION_STYLE = ':'
# MODIFIED Arrow Style
ARROW_STYLE = dict(color=ELECTRON_COLOR, ha='center', va='center', fontsize=12, fontweight='bold') # Was 10, normal
SUPERSCRIPT_MAP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
UP_ARROW = '↑'
DOWN_ARROW = '↓'

def plot_mo_diagram(
    aos_A: List[AtomicOrbital],
    aos_B: List[AtomicOrbital],
    mos: List[MolecularOrbital], # Filled MOs, correctly ordered
    atom_A_label: str,
    atom_B_label: str,
    title: str,
    bond_order: float,
    magnetic_behavior: str,
    output_filename: str = "mo_diagram.svg"
):
    """Generates and saves an MO diagram using Matplotlib."""

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.set_xticks([])
    ax.tick_params(axis='y', direction='in', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['left'].set_color('darkgrey')
    ax.set_ylabel("Energy (eV)", fontsize=10, labelpad=10) #

    # --- Plotting Positions ---
    ao_A_x = 0.15
    ao_B_x = 0.85
    mo_x_center = 0.50
    level_width = 0.08
    ao_level_width = level_width * 1.2
    # MODIFIED Spacing
    degenerate_slot_spacing = level_width * 1.5 # Was 1.1
    electron_vertical_offset_factor = 0.015
    # MODIFIED Electron Offset
    electron_x_offset = level_width * 0.35 # Was 0.25

    # --- Determine Energy Range ---
    all_ao_energies = [ao.energy for ao in aos_A + aos_B]
    all_mo_energies = [mo.energy for mo in mos]
    if not all_ao_energies and not all_mo_energies:
        min_energy, max_energy = -20, 10
    else:
        all_energies = all_ao_energies + all_mo_energies
        min_e = min(all_energies) if all_energies else -10
        max_e = max(all_energies) if all_energies else 10
        padding = max(2.0, abs(max_e - min_e) * 0.12)
        min_energy = min_e - padding
        max_energy = max_e + padding * 1.8

    ax.set_ylim(min_energy, max_energy)
    energy_range = max_energy - min_energy
    if energy_range == 0: energy_range = 1.0
    electron_vertical_offset_abs = electron_vertical_offset_factor * energy_range #

    # --- Helper to Plot Orbitals ---
    def plot_orbital_level(x_center, energy, level_width, color, linewidth, label=None, label_pos='right', va='center'):
        line_start = x_center - level_width / 2
        line_end = x_center + level_width / 2
        ax.hlines(energy, line_start, line_end, color=color, lw=linewidth)
        if label:
            text_x = line_end + 0.018 if label_pos == 'right' else line_start - 0.018
            ha = 'left' if label_pos == 'right' else 'right'
            ax.text(text_x, energy, label, ha=ha, va=va, color=color, fontsize=8.5) #

    # --- Helper to Plot Electrons ---
    def plot_electrons(x_center, y_base, electrons):
        if not electrons: return
        y_pos = y_base + electron_vertical_offset_abs
        if len(electrons) == 1:
            arrow = UP_ARROW
            ax.text(x_center, y_pos, arrow, **ARROW_STYLE) # Uses updated ARROW_STYLE
        elif len(electrons) == 2:
            ax.text(x_center - electron_x_offset, y_pos, UP_ARROW, **ARROW_STYLE) # Uses updated offset & style
            ax.text(x_center + electron_x_offset, y_pos, DOWN_ARROW, **ARROW_STYLE) # Uses updated offset & style

    # --- Plot AOs ---
    # (No changes needed in this specific section based on request)
    ao_levels_A: Dict[float, List[AtomicOrbital]] = {}
    for ao in aos_A:
        key = round(ao.energy, 4)
        if key not in ao_levels_A: ao_levels_A[key] = []
        ao_levels_A[key].append(ao)
    for energy, aos_at_level in sorted(ao_levels_A.items()):
        symbols = "/".join(sorted(list(set(ao.symbol for ao in aos_at_level))))
        plot_orbital_level(ao_A_x, energy, ao_level_width, AO_COLOR, AO_LINE_WIDTH, label=symbols, label_pos='left') #

    ao_levels_B: Dict[float, List[AtomicOrbital]] = {}
    for ao in aos_B:
        key = round(ao.energy, 4)
        if key not in ao_levels_B: ao_levels_B[key] = []
        ao_levels_B[key].append(ao)
    for energy, aos_at_level in sorted(ao_levels_B.items()):
        symbols = "/".join(sorted(list(set(ao.symbol for ao in aos_at_level))))
        plot_orbital_level(ao_B_x, energy, ao_level_width, AO_COLOR, AO_LINE_WIDTH, label=symbols, label_pos='right') #

    # --- Plot MOs ---
    mo_levels: Dict[float, List[MolecularOrbital]] = {}
    for mo in mos:
        key = round(mo.energy, 4)
        if key not in mo_levels: mo_levels[key] = []
        mo_levels[key].append(mo)

    plotted_connections = set()
    processed_groups_for_labels: Set[Tuple[float, str]] = set() # (energy_key, group_id_or_symbol)

    for energy_key, mos_at_level in sorted(mo_levels.items()):
        energy = mos_at_level[0].energy # Use actual energy
        mos_at_level.sort(key=lambda m: m.symbol) # Sort by exact symbol (e.g., _1 before _2)
        num_orbitals = len(mos_at_level)
        # Use MODIFIED degenerate_slot_spacing
        total_width = (num_orbitals - 1) * degenerate_slot_spacing if num_orbitals > 1 else 0
        start_x = mo_x_center - total_width / 2

        # Plot slots and electrons
        for slot_index, mo in enumerate(mos_at_level):
            slot_center_x = start_x + slot_index * degenerate_slot_spacing
            mo_line_color = MO_COLOR
            # Use MODIFIED colors
            if mo.is_bonding is True: mo_line_color = BONDING_COLOR
            elif mo.is_bonding is False: mo_line_color = ANTIBONDING_COLOR
            elif mo.is_nonbonding: mo_line_color = NONBONDING_COLOR

            # Plot orbital line (using MODIFIED MO_LINE_WIDTH)
            ax.hlines(energy, slot_center_x - level_width / 2, slot_center_x + level_width / 2, color=mo_line_color, lw=MO_LINE_WIDTH)
            # Plot electrons (uses updated plot_electrons helper)
            plot_electrons(slot_center_x, energy, mo.electrons)

            # Draw Connecting Lines
            # (No changes needed in this specific section based on request)
            for origin_ao in mo.origin_AOs:
                ao_x = ao_A_x if origin_ao.atom_label == 'A' else ao_B_x
                ao_energy = origin_ao.energy
                ao_key = (round(ao_energy, 3), round(ao_x, 2), origin_ao.symbol)
                mo_key = (round(energy, 3), round(slot_center_x, 2), mo.symbol)
                conn_key = tuple(sorted((ao_key, mo_key)))
                if conn_key not in plotted_connections:
                    ax.plot([ao_x, slot_center_x], [ao_energy, energy], ls=CONNECTION_STYLE, color=LINE_COLOR, lw=0.8, zorder=-1)
                    plotted_connections.add(conn_key) #

        # Add MO Level Labels (handle degeneracy)
        # (No specific changes requested for label logic, but it uses updated positions/colors)
        group_representatives: Dict[str, MolecularOrbital] = {} # group_id -> first MO in group
        group_electrons: Dict[str, int] = {}
        group_x_coords: Dict[str, List[float]] = {}

        for slot_index, mo in enumerate(mos_at_level):
            group_id = mo.degeneracy_group if mo.degeneracy_group else mo.symbol # Unique ID
            slot_center_x = start_x + slot_index * degenerate_slot_spacing

            if group_id not in group_representatives:
                group_representatives[group_id] = mo
                group_electrons[group_id] = 0
                group_x_coords[group_id] = []
            group_electrons[group_id] += mo.electron_count
            group_x_coords[group_id].append(slot_center_x) #

        # Now create the labels
        label_x_positions = []
        labels_to_render = []

        for group_id, rep_mo in group_representatives.items():
            label_key = (energy_key, group_id)
            if label_key not in processed_groups_for_labels:
                # Calculate center X for the group label
                coords = group_x_coords[group_id]
                label_x = (min(coords) + max(coords)) / 2 + level_width / 2 + 0.018 # Right of group center
                label_x_positions.append(label_x)

                total_electrons = group_electrons[group_id]
                # Use the representative MO's symbol, but simplify it for display
                base_symbol = rep_mo.symbol
                if rep_mo.degeneracy_group:
                     base_symbol = base_symbol.rsplit('_', 1)[0]

                electron_str = str(total_electrons).translate(SUPERSCRIPT_MAP) if total_electrons > 0 else ""
                full_label = f"{base_symbol}{electron_str}"
                labels_to_render.append({'x': label_x, 'text': full_label})
                processed_groups_for_labels.add(label_key) #

        # Render labels, aligning them to the rightmost position to avoid overlap
        max_label_x = max(label_x_positions) if label_x_positions else mo_x_center + level_width / 2 + 0.018
        for lbl in labels_to_render:
             ax.text(max_label_x, energy, lbl['text'], ha='left', va='center', color=MO_COLOR, fontsize=8.5)


    # --- Add Atom Labels and Title ---
    label_y_pos = max_energy - padding * 0.15
    ax.text(ao_A_x, label_y_pos, f"Atom A\n{atom_A_label}", ha='center', va='top', fontsize=10, fontweight='medium')
    ax.text(ao_B_x, label_y_pos, f"Atom B\n{atom_B_label}", ha='center', va='top', fontsize=10, fontweight='medium')
    ax.text(mo_x_center, label_y_pos, "MOs", ha='center', va='top', fontsize=10, fontweight='medium') #

    bond_order_str = f"{bond_order:.1f}".rstrip('0').rstrip('.') if '.' in f"{bond_order:.1f}" else str(int(bond_order))
    formatted_behavior = magnetic_behavior.replace('e-', 'e⁻')
    full_title = f"MO Diagram for {title}\nBond Order = {bond_order_str} | {formatted_behavior}"
    # MODIFIED Title Fontsize/Padding
    ax.set_title(full_title, fontsize=14, pad=30, fontweight='bold') # Was 12, pad 25

    # --- Final Touches & Save ---
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    try:
        output_format = 'svg'
        if '.' in output_filename:
             ext = output_filename.split('.')[-1].lower()
             if ext in ['png', 'svg', 'pdf', 'jpg', 'jpeg']: output_format = ext
             else: output_filename = f"{output_filename.rsplit('.', 1)[0]}.svg"
        else: output_filename = f"{output_filename}.svg"
        plt.savefig(output_filename, bbox_inches='tight', dpi=150, format=output_format)
        print(f"Diagram saved to '{output_filename}'")
    except Exception as e:
        print(f"Error saving plot to '{output_filename}': {e}", file=sys.stderr)
    finally:
        plt.close(fig) #