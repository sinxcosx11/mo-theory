from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from .atoms import get_ao_energy, SIGMA_PI_INVERSION_Z # Assuming SIGMA_PI_INVERSION_Z might still be useful
import math
import sys # For warnings
import periodictable # Need this for Z lookup if needed

# Constants (Keep previous values or adjust as needed)
DEFAULT_BETA = -3.0 # eV
ENERGY_DIFF_THRESHOLD = 11.0 # eV (User adjusted)

# --- AtomicOrbital Class ---
# (No changes needed here for now, but adding symmetry labels could be a future enhancement)
@dataclass(frozen=True, order=True)
class AtomicOrbital:
    """Represents an atomic orbital with its properties."""
    energy: float = field(compare=True) # Use energy for primary sorting
    n: int = field(compare=False)
    l: int = field(compare=False) # 0=s, 1=p, 2=d, ...
    symbol: str = field(compare=False) # e.g., "2s", "2p", "3d"
    atom_label: str = field(compare=False) # 'A' or 'B'
    element_symbol: str = field(compare=False) # e.g., 'N', 'O', 'Fe'

    def __post_init__(self):
        if self.l not in [0, 1, 2]:
            raise ValueError(f"AtomicOrbital currently only supports l=0(s), l=1(p), l=2(d), got l={self.l} for {self.symbol}")

# --- MolecularOrbital Class ---
# (No changes needed here)
@dataclass
class MolecularOrbital:
    """Represents a molecular orbital resulting from AO interactions."""
    energy: float = field(compare=True)
    symbol: str = field(compare=False)
    is_bonding: Optional[bool] = field(compare=False)
    origin_AOs: Tuple[AtomicOrbital, ...] = field(compare=False)
    capacity: int = field(default=2, compare=False)
    electrons: List[int] = field(default_factory=list, compare=False)
    degeneracy_group: Optional[str] = field(default=None, compare=False)

    @property
    def electron_count(self) -> int: return len(self.electrons)
    @property
    def is_antibonding(self) -> bool: return self.is_bonding is False
    @property
    def is_nonbonding(self) -> bool: return self.is_bonding is None
    def add_electron(self, spin: int):
        if not (spin == 1 or spin == -1): raise ValueError("Spin must be +1 or -1.")
        if self.electron_count >= self.capacity: raise ValueError(f"MO {self.symbol} full.")
        if self.capacity == 2 and spin in self.electrons: raise ValueError(f"Pauli violation in {self.symbol}.")
        self.electrons.append(spin)
        self.electrons.sort(reverse=True)
    def __lt__(self, other):
        if not isinstance(other, MolecularOrbital): return NotImplemented
        if not math.isclose(self.energy, other.energy, rel_tol=1e-6): return self.energy < other.energy
        return self.symbol < other.symbol
    def __repr__(self):
        elec_str = "".join(['↑' if s == 1 else '↓' for s in self.electrons])
        bond_str = "Bond" if self.is_bonding else ("Anti*" if self.is_antibonding else "NonB")
        return f"MO({self.symbol}, E={self.energy:.2f}eV, Type={bond_str}, Occ=[{elec_str}], DegGrp={self.degeneracy_group})"

# --- Helper Functions for Symmetry (Simplified) ---

def get_symmetry_type(ao: AtomicOrbital) -> str:
    """Approximates AO symmetry as 'sigma' or 'pi'."""
    if ao.l == 0: # s orbital
        return 'sigma'
    elif ao.l == 1: # p orbital
        # Approximation: Assume only one p orbital acts as sigma (pz)
        # For grouping purposes here, we might initially label all p as potentially sigma OR pi
        # The interaction logic later will clarify. Let's return 'p' and handle below.
        return 'p'
    elif ao.l == 2: # d orbital
        # d orbitals have more complex symmetries (sigma, pi, delta)
        # Simplification: label as 'd' and handle interactions specifically
        return 'd'
    else:
        return 'unknown'

def get_mo_symmetry_label(ao1: AtomicOrbital, ao2: Optional[AtomicOrbital] = None, interaction_type: str = 'sigma') -> str:
    """Determines g/u label for homonuclear diatomics."""
    if ao2 is None or ao1.element_symbol != ao2.element_symbol:
        return "" # No g/u for heteronuclear or NB

    l1 = ao1.l
    if interaction_type == 'sigma':
        # s+s sigma -> g; pz+pz sigma -> g; dz2+dz2 sigma -> g (assuming correct AO pairing)
         if l1 == 0 or l1 == 1 or l1 == 2: # Crude check
              return "g"
    elif interaction_type == 'pi':
        # px+px pi -> u; py+py pi -> u; dxz+dxz pi -> u etc.
         if l1 == 1 or l1 == 2: # Crude check
              return "u"
    elif interaction_type == 'delta':
         # dxy+dxy delta -> g; dx2y2+dx2y2 delta -> g
         if l1 == 2:
              return "g"
    return ""

def assign_antibonding_symmetry(bonding_sym: str) -> str:
    """Assigns antibonding g/u label."""
    if bonding_sym == "g": return "u"
    if bonding_sym == "u": return "g"
    return ""

# ------
# if it's called only with symmetry-matched AOs (sigma-sigma or pi-pi)

def combine_atomic_orbitals(ao1: AtomicOrbital, ao2: AtomicOrbital, interaction_type: str, beta: float = DEFAULT_BETA) -> List[MolecularOrbital]:
    """
    Calculates MOs for a specific interaction type (sigma or pi)
    between two AOs assumed to have compatible symmetry and energy.

    Args:
        ao1: First AtomicOrbital.
        ao2: Second AtomicOrbital.
        interaction_type: 'sigma' or 'pi' (or 'delta' if implemented).
        beta: Base interaction parameter (negative eV).

    Returns:
        List containing bonding and antibonding MolecularOrbital(s).
    """
    mos = []
    e1, e2 = ao1.energy, ao2.energy
    l1, l2 = ao1.l, ao2.l
    sym1, sym2 = ao1.symbol, ao2.symbol
    is_homonuclear = (ao1.element_symbol == ao2.element_symbol)

    beta_factors = {'sigma': 1.0, 'pi': 0.8, 'delta': 0.6} # Relative strength
    if interaction_type not in beta_factors:
        print(f"Warning: Unknown interaction type '{interaction_type}' in combine_atomic_orbitals.", file=sys.stderr)
        return []

    beta_eff = beta * beta_factors[interaction_type]

    # Optional: Adjust beta_eff slightly if l1 != l2 (e.g., s-p sigma mixing)
    if l1 != l2 and interaction_type == 'sigma':
        beta_eff *= 0.9 # Slightly weaker interaction for mixed sigma

    # Calculate energy splitting
    energy_diff = abs(e1 - e2)
    effective_beta_sq = beta_eff**2
    delta_e_sq = energy_diff**2
    sqrt_term_val = max(0, delta_e_sq + 4 * effective_beta_sq)
    interaction_term = math.sqrt(sqrt_term_val)

    e_bonding = 0.5 * (e1 + e2 - interaction_term)
    e_antibonding = 0.5 * (e1 + e2 + interaction_term)

    # Ensure reasonable energy ordering relative to AOs
    min_ao_e, max_ao_e = min(e1, e2), max(e1, e2)
    e_bonding = min(e_bonding, min_ao_e - 0.01) # Must be below lowest AO
    e_antibonding = max(e_antibonding, max_ao_e + 0.01) # Must be above highest AO

    # --- MO Labeling ---
    # Base symbol (σ or π)
    mo_base_symbol = "σ" if interaction_type == "sigma" else \
                     "π" if interaction_type == "pi" else \
                     "δ" if interaction_type == "delta" else "?"

    # Symmetry label (g/u for homonuclear)
    sym_label_bond = get_mo_symmetry_label(ao1, ao2, interaction_type)
    sym_label_anti = assign_antibonding_symmetry(sym_label_bond)

    # Basis label (e.g., "1s", "2p", or mixed "1s-2p")
    # Use individual symbols if different, common symbol if same
    ao_basis_label = f"{sym1}" if sym1 == sym2 else f"{sym1}-{sym2}"
    # Alternative: Just use the principal quantum number and type if l matches, e.g. "2p"
    if l1 == l2:
         ao_basis_label = f"{ao1.n}{'spdf'[l1]}"

    # Construct final symbols
    final_sym_bond = f"{mo_base_symbol}{sym_label_bond}({ao_basis_label})"
    final_sym_anti = f"{mo_base_symbol}*{sym_label_anti}({ao_basis_label})"

    # Create MOs
    mos.append(MolecularOrbital(symbol=final_sym_bond, energy=e_bonding, is_bonding=True, origin_AOs=(ao1, ao2)))
    mos.append(MolecularOrbital(symbol=final_sym_anti, energy=e_antibonding, is_bonding=False, origin_AOs=(ao1, ao2)))

    return mos



# --- MAJOR REWRITE: `generate_molecular_orbitals` ---

def generate_molecular_orbitals(atom_A_data: Dict, atom_B_data: Dict) -> Tuple[List[AtomicOrbital], List[AtomicOrbital], List[MolecularOrbital]]:
    """
    Generates atomic and molecular orbitals using a symmetry-focused approach.
    Prioritizes like-orbital interactions before mixing.
    """
    all_mos: List[MolecularOrbital] = []
    processed_ao_ids = set() # Keep track of AOs already used in MO formation

    # --- 1. Create AO instances ---
    # (Helper function to create multiple AOs for p/d levels - same as before)
    def create_aos(atom_data, label) -> List[AtomicOrbital]:
        aos_list = []
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        for orb_type in atom_data['valence_orbital_types']:
            energy = get_ao_energy(atom_data['symbol'], orb_type) # Error handling in get_ao_energy
            n = int(orb_type[0])
            l = l_map[orb_type[1]]
            count = 1
            if l == 1: count = 3
            elif l == 2: count = 5
            for i in range(count):
                 unique_symbol = f"{orb_type}_{i}" if count > 1 else orb_type
                 aos_list.append(AtomicOrbital(energy=energy, n=n, l=l, symbol=unique_symbol,
                                               atom_label=label, element_symbol=atom_data['symbol']))
        return aos_list

    aos_A = create_aos(atom_A_data, 'A')
    aos_B = create_aos(atom_B_data, 'B')
    all_original_aos = sorted(aos_A + aos_B, key=lambda ao: ao.energy) # Keep original list for plotting

    print(f"\nDebug: Initial AOs (Count: {len(all_original_aos)})")

    # --- 2. Identify Sigma Interactions (Revised Logic) ---
    print("\nProcessing Sigma interactions...")

    # Identify potential sigma AOs (s and ONE p per atom)
    sigma_aos_A_candidates = sorted([ao for ao in aos_A if ao.l == 0 or ao.l == 1], key=lambda ao: (ao.l, ao.n))
    sigma_aos_B_candidates = sorted([ao for ao in aos_B if ao.l == 0 or ao.l == 1], key=lambda ao: (ao.l, ao.n))
    sigma_aos_A = []
    sigma_p_ao_A_id = None
    for ao in sigma_aos_A_candidates:
        if ao.l == 0: sigma_aos_A.append(ao)
        elif ao.l == 1 and sigma_p_ao_A_id is None:
            sigma_aos_A.append(ao); sigma_p_ao_A_id = id(ao)
    sigma_aos_B = []
    sigma_p_ao_B_id = None
    for ao in sigma_aos_B_candidates:
        if ao.l == 0: sigma_aos_B.append(ao)
        elif ao.l == 1 and sigma_p_ao_B_id is None:
            sigma_aos_B.append(ao); sigma_p_ao_B_id = id(ao)

    print(f"Debug: Selected Sigma AOs: A={[(ao.symbol, ao.energy) for ao in sigma_aos_A]}, B={[(ao.symbol, ao.energy) for ao in sigma_aos_B]}")

    # --- REVISED SIGMA PAIRING ---
    # Prioritize matching n and l first (e.g., 2s-2s, 2p-2p)

    processed_local_sigma_ids = set() # Track AOs used in this sigma phase

    # Try pairing like-orbitals first (s-s, p-p)
    print("  Attempting like-orbital sigma pairing...")
    for ao_a in sigma_aos_A:
        if id(ao_a) in processed_local_sigma_ids: continue
        for ao_b in sigma_aos_B:
            if id(ao_b) in processed_local_sigma_ids: continue
            # Check for same l-value and energy proximity
            if ao_a.l == ao_b.l and abs(ao_a.energy - ao_b.energy) <= ENERGY_DIFF_THRESHOLD:
                print(f"    Combining Like-Sigma: {ao_a.symbol} ({ao_a.energy:.2f}) + {ao_b.symbol} ({ao_b.energy:.2f})")
                new_mos = combine_atomic_orbitals(ao_a, ao_b, 'sigma')
                all_mos.extend(new_mos)
                processed_local_sigma_ids.add(id(ao_a))
                processed_local_sigma_ids.add(id(ao_b))
                break # Move to next ao_a once ao_b is paired

    # Try s-p mixing only between remaining unused sigma AOs
    print("  Attempting s-p sigma mixing for remaining AOs...")
    remaining_sigma_A = [ao for ao in sigma_aos_A if id(ao) not in processed_local_sigma_ids]
    remaining_sigma_B = [ao for ao in sigma_aos_B if id(ao) not in processed_local_sigma_ids]

    potential_mix_pairs = []
    for ao_a in remaining_sigma_A:
        for ao_b in remaining_sigma_B:
             # Ensure mixing (l values different) and energy proximity
             if ao_a.l != ao_b.l and abs(ao_a.energy - ao_b.energy) <= ENERGY_DIFF_THRESHOLD:
                  potential_mix_pairs.append( (abs(ao_a.energy - ao_b.energy), ao_a, ao_b) )

    potential_mix_pairs.sort() # Closest energy difference first

    for diff, ao_a, ao_b in potential_mix_pairs:
        if id(ao_a) not in processed_local_sigma_ids and id(ao_b) not in processed_local_sigma_ids:
            print(f"    Combining Mixed-Sigma: {ao_a.symbol} ({ao_a.energy:.2f}) + {ao_b.symbol} ({ao_b.energy:.2f})")
            new_mos = combine_atomic_orbitals(ao_a, ao_b, 'sigma')
            all_mos.extend(new_mos)
            processed_local_sigma_ids.add(id(ao_a))
            processed_local_sigma_ids.add(id(ao_b))
            # Allow an AO to participate in mixing only once? For simplicity, yes.
            break # Break inner loop once ao_a finds a partner

    # Update main processed list
    processed_ao_ids.update(processed_local_sigma_ids)
    # --- END REVISED SIGMA PAIRING ---


    # --- 3. Identify Pi Interactions & Non-Bonding Pi ---
    # (Keep the Pi logic from the previous version - it seemed okay for CO)
    print("\nProcessing Pi interactions...")
    pi_aos_A = [ao for ao in aos_A if ao.l == 1 and id(ao) != sigma_p_ao_A_id]
    pi_aos_B = [ao for ao in aos_B if ao.l == 1 and id(ao) != sigma_p_ao_B_id]
    print(f"Debug: Potential Pi AOs: A={[(ao.symbol, ao.energy) for ao in pi_aos_A]}, B={[(ao.symbol, ao.energy) for ao in pi_aos_B]}")

    if pi_aos_A and pi_aos_B:
        repr_pi_ao_a = next((ao for ao in pi_aos_A if id(ao) not in processed_ao_ids), None)
        repr_pi_ao_b = next((ao for ao in pi_aos_B if id(ao) not in processed_ao_ids), None)

        if repr_pi_ao_a and repr_pi_ao_b and abs(repr_pi_ao_a.energy - repr_pi_ao_b.energy) <= ENERGY_DIFF_THRESHOLD:
            print(f"  Combining Pi Level based on: {repr_pi_ao_a.symbol} + {repr_pi_ao_b.symbol}")
            pi_mos_set = combine_atomic_orbitals(repr_pi_ao_a, repr_pi_ao_b, 'pi')
            if len(pi_mos_set) == 2:
                bonding_mo, antibonding_mo = pi_mos_set[0], pi_mos_set[1]
                basis_label = bonding_mo.symbol.split('(')[-1].split(')')[0]
                pi_bond_group = f"pi({basis_label})"
                pi_anti_group = f"pi*({basis_label})"
                num_pi_pairs_to_form = 2
                pi_aos_a_available = [ao for ao in pi_aos_A if id(ao) not in processed_ao_ids]
                pi_aos_b_available = [ao for ao in pi_aos_B if id(ao) not in processed_ao_ids]
                for i in range(min(num_pi_pairs_to_form, len(pi_aos_a_available), len(pi_aos_b_available))):
                     b_sym = f"{bonding_mo.symbol.split('(')[0]}_{i+1}({basis_label})"
                     a_sym = f"{antibonding_mo.symbol.split('(')[0]}_{i+1}({basis_label})"
                     all_mos.append(MolecularOrbital(symbol=b_sym, energy=bonding_mo.energy, is_bonding=True, origin_AOs=(repr_pi_ao_a, repr_pi_ao_b), degeneracy_group=pi_bond_group))
                     all_mos.append(MolecularOrbital(symbol=a_sym, energy=antibonding_mo.energy, is_bonding=False, origin_AOs=(repr_pi_ao_a, repr_pi_ao_b), degeneracy_group=pi_anti_group))
                     processed_ao_ids.add(id(pi_aos_a_available[i]))
                     processed_ao_ids.add(id(pi_aos_b_available[i]))

    for ao in pi_aos_A + pi_aos_B: # Non-bonding Pi
         if id(ao) not in processed_ao_ids:
              print(f"  Creating Non-Bonding Pi MO from: {ao.symbol} ({ao.energy:.2f})")
              nb_symbol = f"nb_pi({ao.symbol.split('_')[0]})"
              all_mos.append(MolecularOrbital(symbol=nb_symbol, energy=ao.energy, is_bonding=None, origin_AOs=(ao,)))
              processed_ao_ids.add(id(ao))


    # --- 4. Handle Remaining AOs as Non-Bonding ---
    # (Same as before)
    print("\nProcessing Remaining Non-Bonding AOs...")
    all_initial_aos = aos_A + aos_B
    for ao in all_initial_aos:
        if id(ao) not in processed_ao_ids:
            print(f"  Creating Non-Bonding MO from: {ao.symbol} ({ao.energy:.2f})")
            nb_symbol = f"nb({ao.symbol.split('_')[0]})"
            if ao.l == 2: nb_symbol = f"nb_d({ao.symbol.split('_')[0]})"
            all_mos.append(MolecularOrbital(symbol=nb_symbol, energy=ao.energy, is_bonding=None, origin_AOs=(ao,)))
            processed_ao_ids.add(id(ao))

    # --- 5. Final Sort & Return ---
    # (Same as before)
    all_mos.sort()
    print("Warning: Sigma/Pi inversion logic skipped in this refactor. Review if needed.")
    print(f"\nDebug: Final MOs generated ({len(all_mos)}):")
    for i, mo in enumerate(all_mos):
         origin_str = ", ".join(f"{ao.symbol}({ao.element_symbol})" for ao in mo.origin_AOs)
         print(f"  {i+1}. {mo.symbol:<20} E={mo.energy:<7.2f} eV  Type: {'Bond' if mo.is_bonding else ('Anti*' if mo.is_antibonding else 'NonB'):<5} Origin: [{origin_str}] DegenGrp: {mo.degeneracy_group}")

    return all_original_aos, all_mos