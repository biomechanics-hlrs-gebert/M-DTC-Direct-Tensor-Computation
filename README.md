# Direct Tensor Computation (DTC)
![Architecture](https://img.shields.io/badge/Architecture-x86-green)
![OS](https://img.shields.io/badge/Linux-64Bit-green)
![version](https://img.shields.io/badge/version-1.0.0-green)

An HPC software for the numerical evaluation of cubic control volumes that are defined as macroscopic finite elements. Resulting in (anisotropic) stiffness matrices.

Intended to compute the mechanical properties of human cancellous bone, a heterogenous, anisotropic structure.

This program formerly was known as *struct process*.

Many thanks to Dr.-Ing. Ralf Schneider who developed the program at its core.
## Usage
The program is invoked like other mpi-parallel applications too. The following is an example that needs embedding in a batch-script. Examples of \*.pbs files for Altair PBS Professional are provided in ```./auxiliaries```

```
mpirun ./bin/dtc_v1.0.0_x86_64 -np 32768 <basename>.meta
```
Testing on a local machine with short turnaround times is recommended with up to 6 cores, depending on your system. 

In general, the tool is capable of running on more than 100.000 x86_64 cores with no limit in wall time. Please be aware that computations on large datasets will need computational power of this magnitude.

## Datasets
... are transferred via file exchange and are not pushed into the repository. 

Computed tomography datasets are fed into the DTC process chain via a meta-file-format. Consisting of a basename-nomenclature and various suffixes, the data are given in their raw, binary format.

## Usage:
The program currently only accepts \*.raw files, given with a proper meta-file and basename.
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

The installation of Open MPI, METIS and PETSc is simplified with the install scripts in ```./lib```.

### Optional: Gnu debugging
I highly recommend the MPI-parallel debugging with forge. Nevertheless, local debugging with gdb works out as well.

[Arm DDT Forge](https://www.arm.com/products/development-tools/server-and-hpc/forge/ddt), [gdb](https://www.gnu.org/software/gdb/), [tmpi](https://github.com/Azrael3000/tmpi), [tmux](https://github.com/tmux/tmux/wiki)

## Set up the Environment

To execute the program, the paths of the libraries must be given as environment variables.
Modify your program according to the requirements of your HPC cluster:
```vim ./auxiliaries/system_environments/<system>.sh```

Then source the script.
```source ./environment.source <system>``` 

Both, compilation and execution of the program need this setup. 

## Build
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
