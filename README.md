# Nanoscale Thermal Transport in Al<sub>x</sub>Ga<sub>1−x</sub>N using Machine Learned Interatomic Potential

This repository contains computational tools and workflows for studying thermal transport properties in Al<sub>x</sub>Ga<sub>1−x</sub>N alloys using the Atomic Cluster Expansion (ACE) machine-learned interatomic potential.

## Overview

This project implements a comprehensive computational framework for:
- **Force constant calculations** using ACE potential at zero temperature
- **Thermal conductivity calculations** via Green-Kubo molecular dynamics (MD)
- **Force constant mapping** using Virtual Crystal Approximation (VCA) for alloy systems
- **Format conversion** between different force constant representations (ShengBTE, TDEP)

## Repository Structure

```
ace_mlip_algan/
├── data/                                   
│   ├── algan_ace_potential.yaml            # ACE potential file
│   └── all_algan_aimd_dft.gzip             # AIMD DFT training data
│
├── zero_kelvin_lattice_dynamics/           # Zero-temperature lattice dynamics
│   ├── calculate_force_constants.py        # Calculate FC2 and FC3 using ACE with Phono3py
│   └── vca_force_constants.ipynb           # Virtual Crystal Approximation for force constants
│
├── force_constant_mapping/                  # Force constant mapping and conversion
│   ├── map_force_constants_2nd_order.py     # Map 2nd-order force constants (FC2)
│   ├── map_force_constants_3rd_order.py     # Map 3rd-order force constants (FC3)
│   ├── utils_fc2.py                         
│   ├── utils_fc3.py                         
│   ├── vca_force_constants_2nd_order.ipynb   # VCA interpolation for FC2
│   ├── vca_force_constants_3rd_order.ipynb   # VCA interpolation for FC3
│   ├── vca_unit_cell.ipynb                 
│   ├── convert_shengbte_to_tdep_fc2.ipynb  
│   ├── convert_tdep_to_shengbte_fc3.ipynb  
│   └── count_nonzero_force_constants.ipynb  
│
└── green_kubo_md/                          # Green-Kubo MD calculations
    ├── calculate_thermal_conductivity.py   # Post-process heat flux data and compute thermal conductivity
    ├── lammps_green_kubo_script.in         # LAMMPS input script for GK-MD
    └── algan_ace_potential.yace            # ACE potential file 
```

## Key Components

### 1. Zero-Temperature Lattice Dynamics (`zero_kelvin_lattice_dynamics/`)

Computes harmonic (FC2) and anharmonic (FC3) force constants using the ACE potential:
- Uses `Phono3py` for generating displacement supercells
- Calculates forces using `PyACECalculator` (ASE interface)
- Outputs force constants in HDF5 format compatible with ShengBTE

**Main script:** `calculate_force_constants.py`

### 2. Force Constant Mapping (`force_constant_mapping/`)

Tools for mapping and averaging force constants in alloy systems:
- **VCA (Virtual Crystal Approximation)**: Combines force constants from AlN and GaN to create Al<sub>x</sub>Ga<sub>1−x</sub>N force constants
- **Format conversion**: Converts between ShengBTE and TDEP force constant formats
- **Analysis**: Counts and analyzes non-zero force constants

**Key scripts:**
- `map_force_constants_2nd_order.py`: Maps FC2 for different Al/Ga compositions
- `map_force_constants_3rd_order.py`: Maps FC3 for different Al/Ga compositions

### 3. Green-Kubo MD (`green_kubo_md/`)

Calculates thermal conductivity using equilibrium MD:
- **LAMMPS simulation**: Runs NVE MD with heat flux calculation
- **Post-processing**: Computes heat current autocorrelation function (HCACF) and thermal conductivity
- **Statistical analysis**: Handles multiple random seeds for error estimation

**Key files:**
- `lammps_green_kubo_script.in`: LAMMPS input for Green-Kubo MD
- `calculate_thermal_conductivity.py`: Post-processes heat flux data and computes thermal conductivity

## Dependencies

### Packages
- `numpy`
- `matplotlib`
- `ase` (Atomic Simulation Environment)
- `phonopy`
- `phono3py`
- `pyace` (ACE potential interface)


### External Software
- **LAMMPS**: For molecular dynamics simulations (with PACE potential support)
- **TDEP**: For temperature-dependent effective potential and phonon calculations at finite temperature

## Usage

### Mapping Force Constants for Alloys

```bash
cd force_constant_mapping
python map_force_constants_2nd_order.py  # For TDEP FC2
python map_force_constants_3rd_order.py  # For TDEP FC3
```

Modify the `fraction_al` variable in the scripts to set the Al composition.

## Data Files

- **ACE Potential**: The trained ACE potential is provided in both YAML (`data/algan_ace_potential.yaml`) and YACE (`green_kubo_md/algan_ace_potential.yace`) formats
- **Training Data**: AIMD DFT data used for training is stored in `data/all_algan_aimd_dft.gzip`

## Citation

## License

## Contact
