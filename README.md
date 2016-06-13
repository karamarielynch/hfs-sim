# Hyperfine Spectrum Simulator

![alt text](https://img.shields.io/badge/License-MIT-blue.svg 'License')
![alt text](https://img.shields.io/badge/Python-2.7-green.svg 'Python version')
![alt text](https://img.shields.io/badge/Tested_on-Mac/Windows-green.svg 'Supported platform')
![alt text](https://img.shields.io/badge/Not_tested_on-Linux-red.svg 'Unsupported platform')

## Purpose
The Hyperfine Spectrum Simulator is a GUI for viewing the atomic hyperfine structure of isotopes. The ability of modify the nuclear spin, hyperfine-structure parameters, linewidth and energy of the atom allows the user to see the effect on the hyperfine structure.  

## Contributors
- Kara Marie Lynch (kara.marie.lynch@cern.ch)

## Features
- Written in Python 2.7 and PyQt4
- Tested on Mac OS X and Windows 8 

## Installation
- Install HFS Simulator in directory of choice

## Dependencies
- [NumPy](http://www.numpy.org/)
- [PyQt4](https://pypi.python.org/pypi/PyQt4)
- [Matplotlib](http://matplotlib.org/)

## Usage
Create `Isotopes.txt` file that contians information on the atomic hyperfine structure that you want to simulate.
Example `Isotopes_Fr.txt` file included for use.
### Basic usage:
`$ python hfs_sim.py Isotopes_Fr.txt`

## Functions:
Simulates hyperfine structure of chosen isotope
- Top plot: Hyperfine spectrum in rest frame of atom (units: Hz)
- Bottom plot: Hyperfine spectrum in laboratory frame (units: cm^-1)
- Outputs: Wavenumber of hyperfine structure peaks (cm^-1) and relative intensities

Options to:
- Choose nuclear spin (I), angular momentum of lower (Jl) and upper (Ju) electronic orbitals
- Choose hyperfine structure A and B parameters for lower and upper states
- Choose linewidth (FWHM) of the spectrum
- Choose to limit hyperfine structure paramters to Au/Al ratio
- Choose voltage (energy, V) of the atom
- Choose wavenumber range and offset to plot
- Click 'Reset range' button to recalculate wavenumber range
