# mo_diagram_generator/cli.py
import argparse
import os

def setup_arg_parser():
    """Sets up the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate Molecular Orbital diagrams for diatomic species.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Show defaults in help
        # MODIFIED: Add epilog with example usage [cite: 1]
        epilog="Example: python -m mo_diagram_generator CN+ --output cn_mo.png"
    )
    parser.add_argument(
        "molecule",
        type=str,
        help="The diatomic molecule or ion to analyze (e.g., 'O2', 'N2-', 'CO', 'LiH', 'CN+'). "
             "For single atom symbols like 'N-', it assumes 'N2-'." # [cite: 1]
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None, # Default filename generated dynamically in main.py
        help="Output filename for the diagram (e.g., 'my_diagram.png'). "
             "Supported formats: png, svg, pdf, jpg. If not specified or invalid, "
             "defaults to '<molecule_name>_MO_diagram.svg'." # [cite: 1]
    )
    # Potential future arguments:
    # parser.add_argument(
    #     "--beta",
    #     type=float,
    #     # default=orbitals.DEFAULT_BETA, # Needs import or passing value
    #     help="Interaction parameter (beta) in eV, controlling MO splitting."
    # )
    # parser.add_argument(
    #     "--energy-threshold",
    #     type=float,
    #     # default=orbitals.ENERGY_DIFF_THRESHOLD, # Needs import or passing value
    #     help="Energy difference threshold (eV) for AO interaction."
    # )

    return parser # [cite: 1]