#------------------------------------------------------------------------------
#> morphometric_evaluation - topology aware DTC scheduler
#
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#> Date:    11.03.2023
#> LastMod: 11.03.2023
#
#> @brief:
#> Schedule the DTC computation based on the morphometric measures of a dataset
#> Done in Python for better flexibility in terms of PBS batch scheduling etc.
#> Python is way more flexible for research scripting, e.g. on the cluster.
#------------------------------------------------------------------------------
import os, struct, argparse, sys, math
import numpy as np
sys.path.insert(1, '/home/geb/00_bone_eval_chain/P_MOD_Python')
#
from mod_decomposition import *
from mod_mechanical import *
from mod_meta_parser import *
#
factors = {}
factors.update({"AVG_NDS_VOX_HEX08_06": 1.3908456310310637})
factors.update({"AVG_NDS_VOX_HEX08_12": 1.2743639355748966})
factors.update({"AVG_NDS_VOX_HEX08_24": 1.2596287755834055})
factors.update({"AVG_NDS_VOX_HEX08_48": 1.2179251259194124})
factors.update({"AVG_NDS_VOX_HEX08_72": 1.1793611129414954})
#  factors.update({"AVG_NDS_VOX_HEX08_96": 0.5021699727716449}) not plausible!
#
factors.update({"FE_NODES_PER_PART_06": 3000})
factors.update({"FE_NODES_PER_PART_12": 6000})
factors.update({"FE_NODES_PER_PART_24": 6000})
factors.update({"FE_NODES_PER_PART_48": 6000})
factors.update({"FE_NODES_PER_PART_72": 6000})
#
# Solid FEs:
DOF = 3
#
# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
def ik8_to_list(file):
    f = open(file, 'rb')
    lst = []
    while True: 
        chunk = f.read(8)
        if not chunk: break 

        lst.append(struct.unpack('q', chunk)[0])

    f.close()

    return(lst)
# -----------------------------------------------------------------------------
#
FMT_STRING="-------------------------------------------------------------------------------"
#
# -----------------------------------------------------------------------------
# User output
# -----------------------------------------------------------------------------
print(
"""
-------------------------------------------------------------------------------
-- Morphometric Evaluation for the Direct Tensor Computation at HLRS       
-------------------------------------------------------------------------------""")

# -----------------------------------------------------------------------------
# Set the parser for the command line arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Just give the meta file.')
parser.add_argument('-metaF', action="store", dest="meta_file", help="Input meta file).")

# -----------------------------------------------------------------------------
# Parse the command line arguments
# -----------------------------------------------------------------------------
cmd_args = parser.parse_args()
 
# -----------------------------------------------------------------------------
# Check if the cmd line args are valid/make sense
# Check for the file containing the necessary information
# -----------------------------------------------------------------------------
if not bool(cmd_args.meta_file):
    print("Please specify an input file.")
    exit(1)

if not os.path.exists(cmd_args.meta_file):
    print("Input file »" + cmd_args.meta_file + "« does not exist.")
    exit(2)

basename = os.path.splitext(cmd_args.meta_file)[0]

dmn_no_file= basename + ".status_preprocess"
vox_file = basename + ".vox"

file_missing = False
if not os.path.exists(vox_file):
    print("Input file »" + vox_file + "« does not exist.")
    file_missing = True

if not os.path.exists(dmn_no_file):
    print("Input file »" + dmn_no_file + "« does not exist.")
    file_missing = True

if file_missing: 
    exit(2)


# # Sample list of integers
# int_list = [42, 13, -7, 0, 66545535]
# f = open('integers.bin', 'wb')
# for ii in int_list:
#     f.write((ii).to_bytes(8, byteorder='little', signed=True))
# f.close()


# Domain numbers (Voxels per domain)

# -----------------------------------------------------------------------------
# Read the information per domain
# -----------------------------------------------------------------------------
vox_list = ik8_to_list(vox_file)
dmn_no_list = ik8_to_list(dmn_no_file)

print("-- Entries in vox_list:     ", len(vox_list))
print("-- Entries in dmn_no_list:  ", len(dmn_no_list))
print(FMT_STRING)

# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")

mi_el_type  =  keywords.get("MICRO_ELMNT_TYPE")
ma_el_order =  keywords.get("MACRO_ELMNT_ORDER")
size_dmn    =  keywords.get("SIZE_DOMAIN")

# -----------------------------------------------------------------------------
# Define time factors (tf_) to adjust the expected runtime
# -----------------------------------------------------------------------------
if mi_el_type[0] == "HEX08":
    no_nds = 8
    tf_nds = 1
elif mi_el_type[0] == "HEX20":
    no_nds = 20
    tf_nds = 20.0/8.0

# 24 dof with ma_el_order == 1 --> 24 iterations with the solver
# 60 dof with ma_el_order == 2 --> 60 iterations with the solver
if int(ma_el_order[0]) == 1:
    no_lc = 24
    tf_ma_el = 1
elif int(ma_el_order[0]) == 2:
    no_lc = 60
    tf_ma_el = 60.0/24.0

print("-- Time factor - nodes:     ", tf_nds)
print("-- Time factor - load cases:", tf_ma_el)
print(FMT_STRING)

# -----------------------------------------------------------------------------
# Get the expected runtimes per domain
# -----------------------------------------------------------------------------
list_of_runtimes = []
#
# Derived with analyze_runtime_to_BVTV
dmn_size_avg=sum([float(i) for i in size_dmn])/3
#
if dmn_size_avg == 0.6:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_06")
    FE_nodes_part = factors.get("FE_NODES_PER_PART_06")
elif dmn_size_avg == 1.2:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_12")
    FE_nodes_part = factors.get("FE_NODES_PER_PART_12")
elif dmn_size_avg == 2.4:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_24")
    FE_nodes_part = factors.get("FE_NODES_PER_PART_24")
elif dmn_size_avg == 4.8:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_48")
    FE_nodes_part = factors.get("FE_NODES_PER_PART_48")
elif dmn_size_avg == 7.2:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_72")
    FE_nodes_part = factors.get("FE_NODES_PER_PART_72")



for dmn in range(len(dmn_no_list)):
    FE_nodes_dmn = vox_list[dmn] * node_factor

    parts_per_domain = np.floor(FE_nodes_dmn/FE_nodes_part)
    
    expected_runtime = 

    total_factors = tf_nds * tf_ma_el 

    list_of_runtimes.append(expected_runtime * total_factors)



print(total_factors)
