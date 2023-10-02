# Casino2Castep

This program reads an electron density (hereafter, simply density) from CASINO stored
as the Fourier components at a set of G-vectors in reciprocal space
and then performs an FFT to obtain the real space density on a grid that can then be read in CASTEP.

## Installation
### Pre-requisites
In order to compile this program, you will require:
- Fortran compiler, currently only the GNU compiler `gfortran` and Intel `ifort` have been tested/supported.
  If you wish to use another compiler, see below.

- `FFTW3` must be installed. In particular, make sure you have installed the double precision version of FFTW3
   and **not** the single (`float`) and long-double precision versions.

- GNU Make (`gmake`). In principle any version of Make should work but compilation has been tested on version 4.4.0.

### Compilation
The easiest way to install is to use the `Makefile`.

1. Clone or download the repository.

2. Open the `Makefile`

3. Set the `F90` variable to your compiler (either `gfortran` or `ifort`).
   Then set the `BUILD` variable to `fast`.

4. Now run the command:
   ```make```

The executable will be located in the root file directory of the repository. If you wish to recompile from scratch, first run
```make clean```
then repeat the above steps.

## Usage
The program is capable of understanding CASTEP input files and parsing them.
Therefore, the same input parameters you use for CASTEP should work here.

The most basic usage is
```castep2casino <casino_file> <castep_seedname>`

`casino_file` contains the CASINO denisty while `castep_seedname` is the seedname of the CASTEP file, that is to say without the file extensions.
At a bare minimum, you will require the cell file.

In order to convert the G-vectors into multiples of the primitive, reciprocal lattice vectors, we must specify the primitive real lattice vectors.
You can either:
- specify them directly in the `lattice_cart` block
- if you are constructing a supercell of the _primitive_ cell, you can use the `supercell_size` keyword in the cell file. For instance, for a 1 x 2 x 3 supercell,
  `supercell_size 1 2 3`.
  Comment this out using an exclamation mark `!`, if you are using the cell for a CASTEP run.

The second thing you must specify is the required CASTEP grid size. When you execute the program, a prompt will come up requiring you to enter the required CASTEP grid size.
You can either do this on a single line with spaces, e.g.
`1 2 3` or enter the first number, then press `Enter`, the second number and so on.

## Licence
Copyright (C) 2023 V Ravindran, S J Clark.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
