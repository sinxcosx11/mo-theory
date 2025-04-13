# mo_diagram_generator/electron_filler.py
from .orbitals import MolecularOrbital
from typing import List, Dict, Set
import sys
import math # For isclose

def fill_electrons(mos: List[MolecularOrbital], total_valence_electrons: int) -> str:
    """
    Fills electrons into MOs according to Aufbau, Pauli, and Hund's rules.
    Modifies the MO objects in the input list directly.
    Handles degenerate orbitals correctly.
    Returns the final electron configuration string.

    Args:
        mos: List of MolecularOrbital objects, ASSUMED SORTED BY ENERGY (with potential sigma/pi swap applied).
        total_valence_electrons: The total number of valence electrons to fill.

    Returns:
        A string representing the final electron configuration (e.g., "σ(2s)² σ*(2s)² π(2p)⁴ σ(2p)²").
        Uses standard MO notation, omits g/u for heteronuclear automatically via MO symbols.
    """
    electrons_remaining = total_valence_electrons
    if electrons_remaining < 0:
        print(f"Warning: Total valence electrons is negative ({total_valence_electrons}). Cannot fill.", file=sys.stderr)
        return "Error: Negative valence electrons"

    # --- Safety Check: Clear any existing electrons ---
    for mo in mos:
        mo.electrons.clear() #

    # --- Group MOs by energy level for applying Hund's rule ---
    filled_mos: List[MolecularOrbital] = [] # Keep track of MOs that received electrons

    # --- Fill electrons level by level (Aufbau Principle) ---
    processed_indices = set()
    i = 0
    while electrons_remaining > 0 and i < len(mos):
        if i in processed_indices:
            i += 1
            continue

        # Identify all MOs at the current energy level (degenerate or accidentally same)
        current_energy = mos[i].energy
        current_level_mos_indices = [idx for idx, mo in enumerate(mos) if math.isclose(mo.energy, current_energy, rel_tol=1e-5, abs_tol=1e-5)] #
        current_level_mos = [mos[idx] for idx in current_level_mos_indices]

        # Mark these indices as processed for the outer loop
        processed_indices.update(current_level_mos_indices)

        # Sort within the level primarily by degeneracy group then symbol for consistent filling order
        current_level_mos.sort(key=lambda m: (m.degeneracy_group or m.symbol, m.symbol))

        # --- Apply Hund's Rule within the current degenerate level ---
        # MODIFIED Loop Logic as per request

        # 1. Fill singly with spin up (+1) until level is half-filled or electrons run out
        electrons_added_stage1 = 0
        for mo_in_level in current_level_mos:
            if electrons_remaining <= 0: break
            # Add spin UP only if the MO is empty and capacity allows
            if mo_in_level.electron_count == 0 and mo_in_level.capacity > 0:
                 try:
                     mo_in_level.add_electron(1) # Add spin up
                     electrons_remaining -= 1
                     electrons_added_stage1 += 1
                     if mo_in_level not in filled_mos: filled_mos.append(mo_in_level)
                 except ValueError as e:
                      print(f"Warning: Error adding initial spin +1 to {mo_in_level.symbol}: {e}", file=sys.stderr)

        if electrons_remaining <= 0: break

        # 2. Pair up with spin down (-1) until level is full or electrons run out
        electrons_added_stage2 = 0
        for mo_in_level in current_level_mos:
             if electrons_remaining <= 0: break
             # Add spin DOWN only if the MO has exactly one electron and capacity allows pairing
             if mo_in_level.electron_count == 1 and mo_in_level.capacity >= 2:
                 try:
                     mo_in_level.add_electron(-1) # Add spin down
                     electrons_remaining -= 1
                     electrons_added_stage2 += 1
                     # MO should already be in filled_mos from stage 1
                 except ValueError as e:
                     print(f"Warning: Error pairing spin -1 in {mo_in_level.symbol}: {e}", file=sys.stderr)

        # Move to the next energy level implicitly via the outer loop incrementing 'i'
        # (The check `if i in processed_indices:` handles skipping already processed levels)

    # --- Final Check ---
    if electrons_remaining > 0:
        print(f"Warning: {electrons_remaining} electrons could not be placed. Total MO capacity might be insufficient or MO generation flawed.", file=sys.stderr)
    elif electrons_remaining < 0:
        print(f"Error: Electron filling resulted in {electrons_remaining} electrons (overfilled). Logic error.", file=sys.stderr)


    # --- Build Configuration String ---
    # (Keep existing logic for building the configuration string)
    config_parts = []
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    processed_groups_at_energy: Dict[float, Set[str]] = {} # Track {rounded_energy: {group_name_or_symbol, ...}}

    for mo in mos: # Iterate through original sorted list
        if mo.electron_count == 0:
            continue # Skip empty MOs

        energy_key = round(mo.energy, 6)
        group_name = mo.degeneracy_group
        effective_group_id = group_name if group_name else mo.symbol

        if energy_key in processed_groups_at_energy and \
           effective_group_id in processed_groups_at_energy[energy_key]:
            continue

        if group_name:
            group_mos = [m for m in mos if m.degeneracy_group == group_name and math.isclose(m.energy, mo.energy, rel_tol=1e-5)]
            group_electrons = sum(m.electron_count for m in group_mos)

            if group_electrons > 0:
                rep_symbol = group_mos[0].symbol
                base_symbol = rep_symbol.rsplit('_', 1)[0] if rep_symbol.endswith(('_1', '_2', '_3')) else rep_symbol
                electron_str = str(group_electrons).translate(superscript_map)
                config_parts.append(f"{base_symbol}{electron_str}")

            if energy_key not in processed_groups_at_energy:
                 processed_groups_at_energy[energy_key] = set()
            processed_groups_at_energy[energy_key].add(effective_group_id)

        else:
            electron_str = str(mo.electron_count).translate(superscript_map)
            config_parts.append(f"{mo.symbol}{electron_str}")
            if energy_key not in processed_groups_at_energy:
                 processed_groups_at_energy[energy_key] = set()
            processed_groups_at_energy[energy_key].add(mo.symbol) #


    return " ".join(config_parts) #