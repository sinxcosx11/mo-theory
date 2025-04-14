# mo_diagram_generator/main.py

# (Keep imports as they are)
import sys
import os
import traceback

try:
    from .cli import setup_arg_parser
    from .utils import parse_molecule_input
    from .atoms import get_atom_data, get_total_valence_electrons
    from .orbitals import generate_molecular_orbitals
    from .electron_filler import fill_electrons
    from .properties import calculate_bond_order, determine_magnetic_behavior
    from .plotter import plot_mo_diagram
except ImportError as e:
     if "__main__" == __name__ and not __package__:
          script_path = os.path.abspath(__file__)
          parent_dir = os.path.dirname(os.path.dirname(script_path))
          module_name = os.path.basename(os.path.dirname(script_path))
          print(f"Import Error: {e}.\n\nIt seems you are running '{os.path.basename(script_path)}' directly.\n"
                f"Navigate one directory up (to '{parent_dir}') and use:\n"
                f"  python -m {module_name} [MOLECULE] [OPTIONS]\n\n"
                f"Example: python -m {module_name} O2\n", file=sys.stderr)
     else:
          print(f"Import Error: {e}. Check package structure and dependencies.", file=sys.stderr)
     sys.exit(1)


def run():
    """Main execution function for the MO Diagram Generator."""
    parser = setup_arg_parser()
    args = parser.parse_args()

    try:
        # 1. Parse Input
        atom1_sym, atom2_sym, charge = parse_molecule_input(args.molecule)
        is_homonuclear = (atom1_sym == atom2_sym)
        charge_str = ""
        if charge == 1: charge_str = "+"
        elif charge == -1: charge_str = "-"
        elif charge > 1: charge_str = f"{charge}+"
        elif charge < -1: charge_str = f"{-charge}-"
        molecule_label_for_title = f"{atom1_sym}{atom2_sym}{charge_str}"
        original_input_label = args.molecule

        print(f"--- Analyzing Species: {molecule_label_for_title} (from input '{original_input_label}') ---")
        print(f"   Type: {'Homonuclear' if is_homonuclear else 'Heteronuclear'} Diatomic")
        print(f"   Charge: {charge:+}")

        # 2. Get Atomic Data
        print("\nFetching atomic data...")
        atom_A_data = get_atom_data(atom1_sym)
        atom_B_data = get_atom_data(atom2_sym)
        print(f"Atom A ({atom1_sym}): Z={atom_A_data['atomic_number']}, Period={atom_A_data['period']}, Group={atom_A_data['group']}, Valence Shells: {atom_A_data['valence_orbital_types']}")
        print(f"Atom B ({atom2_sym}): Z={atom_B_data['atomic_number']}, Period={atom_B_data['period']}, Group={atom_B_data['group']}, Valence Shells: {atom_B_data['valence_orbital_types']}")

        # 3. Calculate Valence Electrons
        total_valence_electrons = get_total_valence_electrons(atom_A_data, atom_B_data, charge)
        print(f"\nTotal Valence Electrons to fill: {total_valence_electrons}")

        # 4. Generate Atomic and Molecular Orbitals
        print("\nGenerating Atomic and Molecular Orbitals...")
        # Adjust unpacking to match the new return signature (2 values instead of 3)
        all_original_aos, mos = generate_molecular_orbitals(atom_A_data, atom_B_data)
        

        if not mos:
             print("\nError: No molecular orbitals were generated.", file=sys.stderr)
             sys.exit(1)

        # 5. Fill Electrons
        print(f"\nFilling {total_valence_electrons} electrons into MOs...")
        electron_config_str = fill_electrons(mos, total_valence_electrons)
        print(f"\nElectron Configuration:\n  {electron_config_str}")

        # 6. Calculate Properties
        bond_order = calculate_bond_order(mos)
        magnetic_behavior = determine_magnetic_behavior(mos)
        print(f"\nCalculated Properties:")
        print(f"  - Bond Order: {bond_order:.2f}")
        print(f"  - Magnetic Behavior: {magnetic_behavior}")

        # 7. Plot the MO Diagram
        output_filename = args.output
        if not output_filename:
            safe_mol_name = "".join(c if c.isalnum() else '_' for c in molecule_label_for_title)
            safe_mol_name = safe_mol_name.replace('+','_plus').replace('-','_minus')
            output_filename = f"{safe_mol_name}_MO_diagram.svg"

        supported_formats = ['.png', '.svg', '.pdf', '.jpg', '.jpeg']
        try:
             file_base, file_ext = os.path.splitext(output_filename)
             if not file_ext or file_ext.lower() not in supported_formats:
                  original_name = output_filename
                  output_filename = f"{file_base}.svg"
                  print(f"\nWarning: Output filename '{original_name}' lacked supported extension. Using '{output_filename}'.")
        except Exception:
             print(f"Warning: Could not parse output filename '{output_filename}'. Using default.", file=sys.stderr)
             safe_mol_name = "".join(c if c.isalnum() else '_' for c in molecule_label_for_title)
             safe_mol_name = safe_mol_name.replace('+','_plus').replace('-','_minus')
             output_filename = f"{safe_mol_name}_MO_diagram.svg"

        
        # Filter the original AOs for plotting based on atom_label
        aos_A_for_plot = [ao for ao in all_original_aos if ao.atom_label == 'A']
        aos_B_for_plot = [ao for ao in all_original_aos if ao.atom_label == 'B']
        

        print(f"\nGenerating plot: '{output_filename}'...")
        plot_mo_diagram(
            # Pass the filtered AO lists
            aos_A=aos_A_for_plot,
            aos_B=aos_B_for_plot,
            mos=mos, # Pass the electron-filled MOs
            atom_A_label=atom1_sym,
            atom_B_label=atom2_sym,
            title=molecule_label_for_title,
            bond_order=bond_order,
            magnetic_behavior=magnetic_behavior,
            output_filename=output_filename
        )

        print("\n--- Analysis Complete ---")

    except ValueError as e:
        print(f"\nInput or Calculation Error: {e}", file=sys.stderr)
        # Print detailed traceback for ValueErrors as well during debugging
        # traceback.print_exc(file=sys.stderr)
        sys.exit(1)
    except ImportError as e:
         print(f"\nDependency Error: {e}. Ensure libraries are installed.", file=sys.stderr)
         sys.exit(1)
    except FileNotFoundError as e:
         print(f"\nFile Error: {e}. Check required data files.", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr) # Print detailed traceback
        print("-" * 60, file=sys.stderr)
        print("Please check input and data.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    run()
