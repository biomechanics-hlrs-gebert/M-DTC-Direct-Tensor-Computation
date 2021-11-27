# Directly Discretizing Tensor Computation (DDTC)
The calculation of structural tensors based on the results of an FEA
analysis with 24 (or 60) appropriate load cases, depending of the order of the macro element.

This program formerly was known as »struct process«

Build the program:    »make«
Create documentation: »make docs«

Many thanks to Dr.-Ing. Ralf Schneider who developed the program at its core.

## Dataformats
Two different approaches of dealing with I/O are used. 
### External data
Computed tomography datasets are fed into the DDTC process chain via a meta-file-format. Consisting of a basename-nomenclature and various suffixes, the data are given in their raw, binary format.

### Internal data
Internally and for general program output, the PureDat format gets used. 