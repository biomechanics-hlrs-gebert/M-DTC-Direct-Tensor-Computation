# Direct Tensor Computation (DTC)
![Architecture](https://img.shields.io/badge/Architecture-x86-green)
![OS](https://img.shields.io/badge/Linux-64Bit-green)
![version](https://img.shields.io/badge/version-1.0.0-green)


The computation of structural tensors based on the results of an FEA analysis with 24 (or 60) appropriate load cases, depending of the order of the macro element.

This program formerly was known as »struct process«

Build the program:    ```make```
Create documentation: ```make docs```

Many thanks to Dr.-Ing. Ralf Schneider who developed the program at its core.
## Usage
For example for testing on julius:
```
mpirun ./bin/dtc_v1.0.0_x86_64 -np 4 <basename>.meta```
```

## Datasets
... are transfered via file exchange and are not pushed into the repository. 

## Dataformats
Two different approaches of dealing with I/O are used. 
### External data
Computed tomography datasets are fed into the DTC process chain via a meta-file-format. Consisting of a basename-nomenclature and various suffixes, the data are given in their raw, binary format.
### Internal data
Internally and for general program output, the PureDat format gets used. 

## Usage:
The program currently only accepts *.vtk files with a proper meta basename.
```./dtc_v1.0.0_x86_64 <basename>.meta```

## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with Bash or Zsh
* GNU Compiler Collection (GCC), especially with gfortran
* An installation of Open-MPI
* Geberts libraries. Managed by: ```./manage_geb-lib.sh```

The program must be compiled with:
* Global integer kind=64Bit, signed
* Meta-format integer kind=64Bit, signed
* MPI integer kind=32Bit
* PETSc index length=32Bit
* METIS index length=64Bit

The installation of Open MPI, METIS and PETSc is simplified with the install scripts of the repository "Overview" of the biomechanics-hlrs-gebert organization @GitHub.

### Optional: Gnu debugging
* [gdb](https://www.gnu.org/software/gdb/)
* [tmpi](https://github.com/Azrael3000/tmpi)
* [tmux](https://github.com/tmux/tmux/wiki)
## Build
It's tested and therefore recommended to build and run the program as follows.
### Set up the Environment
```vim ./auxiliaries/system_environments/<system>.sh```
```source ./environment.source <system>``` 

* Set an architecture/a system
  * Give the absolute base path of your mpi-installation
  * Alternatively give the proper module names of your compute cluster

### Run make:
Build the program:    ```make```
Create documentation: ```make docs```

### Uninstall:
```make clean && rm -r <your program directory>```
## Acknowledgements 
Plain text parsed via [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
## Arbitrary
Use this program at your own risk.