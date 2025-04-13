# Molecular Orbital (MO) Diagram Generator for Diatomic Species

A Python tool to generate qualitative Molecular Orbital diagrams for simple diatomic molecules and ions. It calculates approximate MO energy levels, fills valence electrons, determines bond order and magnetic properties, and generates a visual plot using Matplotlib.

## Overview

This program takes a simple diatomic species input (e.g., "N2", "CO", "HF", "O2-") and performs the following steps:

1.  Parses the input to identify the atoms and charge.
2.  Retrieves valence atomic orbital (AO) data, including approximate energies (VSIPs).
3.  Constructs Molecular Orbitals (MOs) by considering AO symmetry (sigma/pi approximation) and energy proximity. It supports s-p mixing.
4.  Fills the total valence electrons into the generated MOs according to Aufbau, Pauli, and Hund's rules.
5.  Calculates the bond order and determines if the species is paramagnetic or diamagnetic.
6.  Generates and saves a plot of the MO energy level diagram using Matplotlib.

## Features

* Generates MO diagrams for homonuclear and heteronuclear diatomic species.
* Handles simple positive and negative ions (e.g., N2+, O2-).
* Includes basic support for valence s, p, and d orbitals.
* Attempts to model s-p mixing based on energy proximity.
* Calculates and displays Bond Order.
* Determines and displays Magnetic Behavior (Paramagnetic/Diamagnetic).
* Outputs a visual MO diagram plot (SVG default, other formats like PNG, PDF supported).
* Command-line interface for easy usage.

## Requirements

* Python 3.x
* Libraries listed in `requirements.txt`:
    * `matplotlib>=3.0.0`
    * `periodictable>=1.6.0`

## Installation

1.  **Clone the repository:**
    ```bash
    git clone (https://github.com/sinxcosx11/mo-theory/)
    cd mo
    ```
    

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
    *(It's recommended to do this within a Python virtual environment)*

## Usage

Run the script as a Python module from the directory *containing* the `mo` folder 
```bash
python -m mo.main <MOLECULE> [OPTIONS]
