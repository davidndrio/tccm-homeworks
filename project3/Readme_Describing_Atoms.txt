# Molecular Dynamics Project

## Description
This project implements a molecular dynamics simulation in Fortran. It includes:
1. Reading atomic data (coordinates and masses) from an input file.
2. Calculating internuclear distances between all pairs of atoms.

## File Structure
- `read_Natoms.f90`: Function to read the number of atoms from the input file.
- `read_molecule.f90`: Subroutine to read atomic coordinates and masses.
- `compute_distances.f90`: Subroutine to calculate internuclear distances.
- `main_program.f90`: Main program that integrates all functionalities.

## Input File Format
The input file must follow this format:
1. The first line contains the number of atoms (`Natoms`).
2. Each subsequent line contains the x, y, and z coordinates followed by the mass of an atom.
Compile execution:
 gfortran -o main_program read_Natoms.f90 read_molecule.f90 compute_distances.f90 main_program.f90
