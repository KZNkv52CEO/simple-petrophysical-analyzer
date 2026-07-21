# Simple Petrophysical Log Analysis

A lightweight Python tool for basic petrophysical interpretation using LAS or CSV well log data.

## Features
- Reads standard .las (via lasio`) and .csv` files.
- Automatically maps log mnemonics (`GR`, RHOB, `Rt`).
- Computes Density Porosity.
- Calculates shale volume ($V_{sh}$) using a linear Gamma Ray method.
- Performs shale-corrected effective porosity ($\phi_e$) estimation.
- Computes water saturation ($S_w$) using the Archie equation.
- Generates a multi-track well log plot with inverted depth scale.

## Requirements
- numpy
- pandas
- matplotlib
- lasio

## Usage
Run the script from your terminal:
```bash
python petrophysics_analizer.py <path_to_las_or_csv_file>
