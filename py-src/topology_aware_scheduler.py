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
sys.path.insert(1, '/home/geb/00_bone_eval_chain/P_MOD_Python')
#
from mod_decomposition import *
from mod_mechanical import *
from mod_meta_parser import *
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

print("-- Entries in vox_list:   ", len(vox_list))
print("-- Entries in dmn_no_list:", len(dmn_no_list))
print(FMT_STRING)

# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")

for key, value in keywords.items():
    print(key, value)