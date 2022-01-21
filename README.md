# Directly Discretizing Tensor Computation (DDTC)
The calculation of structural tensors based on the results of an FEA
analysis with 24 (or 60) appropriate load cases, depending of the order of the macro element.

This program formerly was known as »struct process«

Build the program:    ```make```
Create documentation: ```make docs```

Many thanks to Dr.-Ing. Ralf Schneider who developed the program at its core.
## Cite the Directly Discretizing Tensor Computation (BibTex):
```
@misc{DDTC,
title={{Directly Discretizing Tensor Computation: A program to compute the apparent stiffness of control volumes of heterogeneous structures.

.}},
howpublished="\url{https://github.com/JoGebert/Directly_Discretizing_Tensor_Computation-DDTC/}"
}
```
## Dataformats
Two different approaches of dealing with I/O are used. 
### External data
Computed tomography datasets are fed into the DDTC process chain via a meta-file-format. Consisting of a basename-nomenclature and various suffixes, the data are given in their raw, binary format.
### Internal data
Internally and for general program output, the PureDat format gets used. 
## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with Bash or Zsh
* GNU Compiler Collection (GCC), especially with gfortran
* An installation of Open-MPI

The program must be compiled with:
* Global integer kind=64Bit, signed
* Meta-format integer kind=64Bit, signed
* MPI integer kind=32Bit
* PETSc index length=64Bit
* METIS index length=64Bit

The installation of Open MPI, METIS and PETSc is simplified with the install scripts of the repository "Overview" of the biomechanics-hlrs-gebert organization @GitHub.
### Optional: Gnu debugging
* [gdb](https://www.gnu.org/software/gdb/)
* [tmpi](https://github.com/Azrael3000/tmpi)
* [tmux](https://github.com/tmux/tmux/wiki)

## External Sources
Plain text headers are parsed via a [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
## Arbitrary
Use this program at your own risk.

