# mo_diagram_generator/properties.py
from .orbitals import MolecularOrbital
from typing import List

def calculate_bond_order(mos: List[MolecularOrbital]) -> float:
    """
    Calculates the bond order from the list of filled MOs.
    Non-bonding electrons do not contribute.
    """
    bonding_electrons = 0
    antibonding_electrons = 0

    for mo in mos:
        if mo.is_bonding is True: # Explicitly bonding
            bonding_electrons += mo.electron_count
        elif mo.is_bonding is False: # Explicitly antibonding
            antibonding_electrons += mo.electron_count
        # Non-bonding (is_bonding is None) are ignored

    bond_order = (bonding_electrons - antibonding_electrons) / 2.0

    # Bond order cannot be negative
    return max(0.0, bond_order)


def determine_magnetic_behavior(mos: List[MolecularOrbital]) -> str:
    """
    Determines if the molecule is paramagnetic (unpaired electrons) or
    diamagnetic (all electrons paired) based on the filled MOs.
    """
    total_unpaired_electrons = 0
    for mo in mos:
        # An MO contributes to paramagnetism if it is singly occupied (has exactly 1 electron).
        if mo.electron_count == 1:
            total_unpaired_electrons += 1
        # Assumes capacity 2 orbitals. If capacity > 2 were implemented,
        # this check would need refinement based on Hund's rule application.

    if total_unpaired_electrons > 0:
        # Use Unicode minus sign for electron charge symbol
        plural = "s" if total_unpaired_electrons > 1 else ""
        return f"Paramagnetic ({total_unpaired_electrons} unpaired e⁻)"
    else:
        return "Diamagnetic"

