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
## Prerequisites
* Working mpi installation with 
  * mpi-compilers
  * mpirun
* Working installation of METIS with 64Bit index length
* Working installation of PETSc with 64Bit index length
### Optional: Gnu debugging
* [gdb](https://www.gnu.org/software/gdb/)
* [tmpi](https://github.com/Azrael3000/tmpi)
* [tmux](https://github.com/tmux/tmux/wiki)

