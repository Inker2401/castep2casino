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

### Known Issues
Currently, compilation will fail if `Make` is run in parallel using the `-j` flag.

## Usage
The program is capable of understanding CASTEP input files and parsing them.
Therefore, the same input parameters you use for CASTEP should work here.

The most basic usage is
```castep2casino <casino_file> [lat_geom_file]```
where `casino_file` contains the CASINO density.

If you do not specify the `lat_geom_file`, the program will look for it in the form of in the form of `seedname.latt_geom`.
The `seedname` will be taken from the `casino_file`, i.e. if the CASINO file is named `QMC.txt`, the program will look for the file `QMC.latt`.

Note that keywords are case-insensitive.

<details><summary>User Options</summary>

### Summary
The `<lat_geom_file>` effectively contains the set of user-defined parameters which are:
- `prim_latt_cart` block - primitive real lattice vectors.
- `castep_grid` = the dimensions of the CASTEP grid.
- `unit_bohr` - unit for real lattice vectors are Bohr. This is an **OPTIONAL** parameter, and otherwise angstroms are used as the input unit.
- `output_file` - The file to use to write the CASTEP output. This is an **OPTIONAL** parameter.

#### Essential
The primitive real lattice vectors are specified as follows:
```
%block prim_latt_cart
a_x a_y a_z
b_x b_y b_z
c_x c_y c_z
%endblock prim_latt_cart
```

The default input unit for the lattice vectors are ANGSTROMS.
You may alternatively put the string `unit_bohr` anywhere in the cell file and then enter the lattice vectors in atomic units (Bohr).

The second thing you must specify is the required CASTEP grid size. This is likewise done in the `lat_geom_file` through the parameter
```
castep_grid <nx> <ny> <nz>
```
For instance, if you needed a grid of size 10 x 20 x 30, you would enter `castep_grid 10 20 30`.

For the output file, if `output_file` is not specified, the program will use the `lat_geom_file` seedname.
For example, if you have a file named `QMC.dat` as the lattice geometry file, then the output file will be `QMC.den_fmt`.
</details>


## Licence
Copyright (C) 2023 V Ravindran

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
