# Partition function calculation for logical failure estimation of toric code in [1]
This repository contains a modified version of `isingZ` which was originally developed by Creighton K. Thomas (creightonthomas@gmail.com) and A. Alan Middleton (aam@syr.edu). The modifications of this repository aim at generating large sample sets of partition function estimates for a uniform disorder distribution and a truncated Gaussian disorder model.
The finite temperature and Ising couplings are chosen to fit the statistical mechanics mapping of stabilizer codes to relate the partition functions to logical error curves as explained in detail in [1].

## Notes on original authors and copyright
```text
// COPYRIGHT NOTICE
// This repository is a modifed version of isingZ, a program to compute partition
//     functions in 2D Ising models
//     for general bond weights and periodic boundary conditions.
// Copyright 2012 by Creighton K. Thomas (creightonthomas@gmail.com),
//                   A. Alan Middleton (aam@syr.edu)
// The development of this software was supported in part by the National Science
// Foundation under grant DMR-1006731.
//
// We have edited the original code by editing the workflow to output txt files for large scale
// experiments, modified the interaction generation and edited the temperature definition in terms of
// the Nishimori temperature.
//
// The program isingZ is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// See http://www.gnu.org/licenses/gpl.txt for more details.
//
```

## Notes on Authors of Modified Version
This repository was created by the authors of [1], Leon Wichette, Hans Hohenfeld, Elie Mounzer and Linnea Grans-Samuelsson, to utilize `isingZ` for the simulation of partition function decoding on the toric code.

Our changes include:
- **Interaction and partition function output to `.txt`**: Enables large scale experimental setup.
- **Disorder sampling**: The generator now supports both a **uniform** distibution and **truncated Gaussian** disorder distribution.
- **Nishimori condition**: The temperature and interaction strengths are determined by the Nishimori condition to enable the statistical mechanical mapping of stabilizer codes.
- **Simplified codebase**: Unused components have been removed to streamline the codebase.
- **Packing of results**: A script called `combine_to_hdf5.py` was added to store results in a single file which acts as input for the post processing pipeline.

For scientific or reproducibility inquiries, please refer to the arXiv paper above or contact the authors listed therein.


## Short directions for compilation and verification of compile
If you have all the libraries available, which might well be the case, simply
  enter 'make'. This will compile codes and place them in the bin directory.
  It will also run short tests.
  If the tests fails, don't worry, yet: the check is for the full precision of the
  result and it is conceivable that different versions of gmp will give very slight
  differences in the last few digits of the long decimal outputs.

Necessary libraries for ising Z include gmp, gmpxx (GNU multiprecision library)..
The disorder realization generator needs gsl (GNU scientific library).
  You can modify the Makefile in src/Z/ and src/generator/ if you need to
  clarify where those libraries sit. (For example, on one of
  the authors' computers, where gmp and gsl are installed in bizarre places,
  the additions to the makefile include:
    INC		= -I/sw/include/gsl -I/sw/include/
    LIBS	= -lgmpxx.4 -lgmp.3 -L/usr/local/lib -L/sw/lib/x86_64/ -lgsl
  )

## Brief description of input format and usage
```text
Conventions for input format for Ising samples:
 See sample files in src/test/ for examples of input files.
 These files include system size and a list of bond weights, two
   bond weights per spin. The weights are indexed by spin coordinates
   and bond directions.
 Diagram of spin layout =

 (   0,   0) (   1,   0) (   2,   0) ... (Lx-1,   0)
 (   0,   1) (   1,   1) (   2,   1) ... (Lx-1,   1)
 (   0,   2) (   1,   2) (   2,   2) ... (Lx-1,   2)
      .           .           .               .
      .           .           .               .
      .           .           .               .
 (   0,Ly-1) (   1,Ly-1) (   2,Ly-1) ... (Lx-1,Ly-1)


 For this spin layout, directions =

      N(0)

  W(3)    E(1)

      S(2)
```

## Usage instructions

The workflow consists of two steps:
1. Generate Ising interactions using `generator_random_bond`
2. Calculate partition functions using `isingZToTxt`
3. Pack results in HDF5 format using `combine_to_hdf5.py`

### Step 1: Generate Ising Interactions

```bash
./bin/generator_random_bond Lx Ly seed probability output_directory [std_deviation]
```

Parameters:
- `Lx`, `Ly`: Lattice dimensions (integers)
- `seed`: Random number generator seed (integer)
- `probability`: Probability parameter for bond generation (double, 0.0 to 0.5)
- `output_directory`: Base directory for output files
- `std_deviation`: (Optional) Standard deviation for Gaussian noise distribution

Example:
```bash
./bin/generator_random_bond 5 5 42 0.1 ./data
```

This generates a 5×5 lattice with probability 0.1 of flipping a bond and saves it to:
`./data/interactionsGaussian/0.100000/5/5/0.000000/42/interaction_lattice.txt`

With truncated Gaussian distributed coupling flipping probability with mean probability 0.1 and standard deviation 0.05:
```bash
./bin/generator_random_bond 5 5 42 0.1 ./data 0.05
```

### Step 2: Calculate Partition Functions

```bash
./bin/isingZToTxt precision Lx Ly seed probability temperature output_directory [std_deviation]
```

Parameters:
- `precision`: Bits of precision for calculations (e.g., 512, 1024, 4096)
- `Lx`, `Ly`: Same lattice dimensions as used in Step 1
- `seed`: Same seed as used in Step 1
- `probability`: Same probability as used in Step 1
- `temperature`: Temperature factor measured in units of Nishimori temperature
- `output_directory`: Same base directory as used in Step 1
- `std_deviation`: (Optional) Same standard deviation as used in Step 1

Example:
```bash
./bin/isingZToTxt 4096 5 5 42 0.1 1.0 ./data
```

This reads the interaction file and outputs the partition functions to:
`./data/resultsGaussian/0.100000/0.000000/5/5/1.000000/4096/42/Z.txt`

With Gaussian noise:
```bash
./bin/isingZToTxt 4096 5 5 42 0.1 1.0 ./data 0.05
```

The output text file contains four tab-separated partition function values that represent four different boundary conditions of the spin lattice:
 - (periodic, periodic)
 - (periodic, anti periodic)
 - (anti periodic, periodic)
 - (anti periodic, anti periodic)


### Step 3: Combine results

To handle the results easier it may be useful for you to pack the generated results in a structured way into a HDF5 file. THis can be achieved by calling:

```bash
python combine_to_hdf5.py results_dir output_file_name
```
Parameters:
- `results_dir`: The directory the results directory was created in (same directory the algorithm scripts were executed with)
- `output_file_name`: The name of the output HDF5 file

Example:
```bash
python combine_to_hdf5.py ./data results.h5
```
This reads the results from `./data/resultsGaussian` and stores the content in `results.h5`.


## Acknowledgments

This work was funded by the German Ministry of Economic Affairs and Climate Action (BMWK) and the
German Aerospace Center (DLR) in project QuDA-KI under grant no. 50RA2206 and through a
Leverhulme-Peierls Fellowship at the University of Oxford, funded by grant no. LIP-2020-014.


## References

[1] Wichette, L., Hohenfeld, H., Mounzer, E., & Grans-Samuelsson, L. (2025). *A partition function
framework for estimating logical error curves in stabilizer codes*.  arXiv preprint
[arXiv:2505.15758](https://arxiv.org/abs/2505.15758).
