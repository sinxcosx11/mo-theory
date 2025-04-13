# mo_diagram_generator/utils.py
import re
from typing import Tuple

def parse_molecule_input(input_str: str) -> Tuple[str, str, int]:
    """
    Parses input like "O2", "N2-", "Fe2+", "CO", "NO+".
    Returns atom symbols (tuple of 2) and charge (int).
    Focuses on diatomic molecules (homo or hetero). Handles simple ions.

    Args:
        input_str: The user-provided string representing the molecule/ion.

    Returns:
        A tuple containing (atom1_symbol, atom2_symbol, charge).

    Raises:
        ValueError: If the input string cannot be parsed into a diatomic species.
    """
    input_str = input_str.strip()
    # Regex to capture:
    # Group 1: First element symbol ([A-Z][a-z]?)
    # Group 2: Optional digit '2' indicating homonuclear diatomic (2)?
    # Group 3: Optional second element symbol ([A-Z][a-z]?) for heteronuclear
    # Group 4: Optional charge sign ([+-])?
    # Group 5: Optional charge magnitude (\d+)? (defaults to 1 if sign exists)
    # Allows spaces before charge, e.g., "O2 2-"
    pattern = r"^([A-Z][a-z]?)(\d)?([A-Z][a-z]?)?\s*([+-])?(\d+)?$"
    match = re.match(pattern, input_str)

    if not match:
        raise ValueError(
            f"Invalid molecule format: '{input_str}'. "
            "Expected formats like 'O2', 'N2-', 'CO', 'LiH', 'CN+'."
        )

    atom1_sym = match.group(1)
    count1 = match.group(2) # Should be '2' for explicit homonuclear like O2
    atom2_sym_explicit = match.group(3) # For heteronuclear like CO
    charge_sign_str = match.group(4)
    charge_mag_str = match.group(5)

    # Determine charge
    charge = 0
    if charge_sign_str:
        magnitude = int(charge_mag_str) if charge_mag_str else 1
        charge = magnitude if charge_sign_str == '+' else -magnitude

    # Determine the two atoms involved
    atom1_final = atom1_sym
    atom2_final = None

    if count1 == '2' and not atom2_sym_explicit:
        # Homonuclear diatomic specified like O2, N2, etc.
        atom2_final = atom1_sym
    elif atom2_sym_explicit and not count1:
        # Heteronuclear diatomic specified like CO, LiH, etc.
        # Check if atom1_sym and atom2_sym_explicit are the same (e.g., "NN" input)
        if atom1_sym == atom2_sym_explicit:
             print(f"Warning: Input '{input_str}' interpreted as homonuclear '{atom1_sym}2'.")
        atom2_final = atom2_sym_explicit
    elif not count1 and not atom2_sym_explicit:
        # Single atom symbol given (e.g., "He", "Fe", "N-").
        # Assume user means the diatomic molecule (He2, Fe2, N2-) for MO diagram.
        print(f"Warning: Input '{input_str}' interpreted as diatomic '{atom1_sym}2' with charge {charge:+}.")
        atom2_final = atom1_sym
    else:
        # Ambiguous or unsupported format (e.g., "O3", "C2H4", "CO2", "O2F")
        raise ValueError(
            f"Cannot interpret '{input_str}' as a simple diatomic molecule. "
            f"Please use formats like 'O2', 'N2-', 'CO'."
        )

    if atom1_final is None or atom2_final is None:
         # Should not happen if logic above is correct, but as a safeguard
         raise ValueError(f"Failed to determine two atoms from input '{input_str}'.")

    # Return the two atom symbols and the charge
    return atom1_final, atom2_final, charge

