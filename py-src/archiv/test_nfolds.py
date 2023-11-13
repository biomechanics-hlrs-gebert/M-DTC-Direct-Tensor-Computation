import sys
import os
import struct
# -----------------------------------------------------------------------------
# Read the binary neighborhood score files
# -----------------------------------------------------------------------------
#
WORKING_DIRECTORY =  os.getcwd()
#
file_path = WORKING_DIRECTORY + "/24er_morph_aware_V1.6.0_HEX08/FH01-1_mu_prod_dtc_topo-morph-HEX08-ME1-24.6fold"

with open(file_path, 'rb') as file:
    file_content = file.read()

    num_floats = len(file_content) // struct.calcsize('d')

    floats = struct.unpack('d' * num_floats, file_content)

    ii = 0
    for f in floats:
        
        if ii == 10:
            break 
            
        print(f)
        ii += 1
        