# DP4

**Copyright (c) 2024 Alexander Howarth, Kristaps Ermanis, Jonathan M. Goodman, Benji Rowlands**

**Distributed under MIT license**

---

## Contents

1. Release Notes
2. Requirements and Setup
3. Usage
4. NMR Description Format
5. Code Organization
6. Help and Contact Details

---

## Release Notes

This repository contains all the code necessary to reproduce the DP4* results in
the paper **"Towards automatically verifying chemical structures: the powerful
combination of <sup>1</sup>H NMR and IR spectroscopy"**.

This package should be used in conjunction with the [Apollo data repository](https://doi.org/10.17863/CAM.110235) in order to reproduce the DP4* results.

Scripts to reproduce the SCC plots from the paper are also available [on GitHub](https://github.com/Goodman-lab/SCC).

Functionality from other releases of DP4 (diastereoisomer generation, automatic
interpretation of raw NMR data, conformational analysis of 5-membered rings,
calculation of DP5 probabilities) has been removed to facilitate easy
reproduction of the results in the paper.

---

## Requirements and Setup

### Requirements

- Python 3.6+
- `scipy` (version 1.14.0)
- `numpy` (version 2.0.1)
- `rdkit` (version 2024.3.4)
- `networkx` (version 3.3)
- `setuptools` (version 72.1.0)

### Setup Instructions

1. **Clone or download the repository**

   Clone or download the repository to your local machine. You can use the following git command to clone the repository:
   ```
   git clone https://github.com/jbr46/DP4.git
   ```
2. **Create and activate a virtual environment**
   Navigate to the directory where you cloned the repository and run the following command:
   ```
   python -m venv .venv
   ```
   This command creates a virtual environment named `.venv` in your DP4 directory.
3. **Activate the virtual environment**
   ```
   source .venv/bin/activate
   ```
4. **Install the package**
   With the virtual environment activated, install DP4 and its dependencies using `pip`:
   ```
   pip install .
   ```
---

## NMR Description Format

**Example NMR file**:

59.58(any),127.88(any),...
4.81(any),7.18(any),...

### Sections:

1. **C shifts** (e.g., 59.58(any)).
2. **H shifts** (e.g., 4.81(any)).

In this work, the shifts are always assumed to be **unassigned**, i.e., each one should be followed by (any).

---

## Code Organization

- **`dp4/`**: The main package directory.
- **`__init__.py`**: Initializes the dp4 package.
- **`ConfPrune.pyx`**: Cython file for conformer alignment and RMSD pruning.
- **`main.py`**: Main script for DP4 workflow.
- **`MacroModel.py`**: Handles MacroModel-specific calculations.
- **`Gaussian.py`**: Handles Gaussian-specific calculations.
- **`NWChem.py`**: Handles NWChem-specific calculations.
- **`DP4.py`**: Python implementation of DP4 probability calculation.
- **`StructureInput.py`**: Cleans structures, generates 3D geometries, reads SMILES and SMARTS formats.
- **`data/`**: Contains data files such as `TMSdata`.

---

## Help and Contact Details

For more information on using the scripts and additional options, refer to the
in-script help by running:
```
pydp4 -h
```
In case of any questions or issues, please email Jonathan Goodman at
jmg11@cam.ac.uk.
