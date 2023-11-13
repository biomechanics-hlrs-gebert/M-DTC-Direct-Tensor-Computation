#------------------------------------------------------------------------------
#> morphometric_evaluation - topology aware DTC scheduler
#
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#> Date:    10.10.2023
#> LastMod: 10.10.2023
#------------------------------------------------------------------------------
import os, struct, argparse, sys, time
import numpy as np
import pandas as pd
from pathlib import Path
import struct
import shutil
#
sys.path.insert(1, '/home/geb/00_bone_eval_chain/P_MOD_Python')
#
from mod_decomposition import *
from mod_mechanical import *
from mod_meta_parser import *
#
#------------------------------------------------------------------------------
# Replace line in a file
#------------------------------------------------------------------------------
def replace_line_in_file(file_path, keyword, new_line):
    # Read the file and store its lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Create a flag to track if a replacement occurred
    replaced = False

    # Iterate through the lines and replace if the keyword is found
    for i in range(len(lines)):
        if keyword in lines[i]:
            lines[i] = new_line + '\n'
            replaced = True

    # If a replacement occurred, write the modified lines back to the file
    if replaced:
        with open(file_path, 'w') as file:
            file.writelines(lines)
    else:
        print(f"Keyword '{keyword}' not found in {file_path}")
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
FMT_SEP="-------------------------------------------------------------------------------"
start_time = time.time()
#
factors = {}
factors.update({"AVG_NDS_VOX_HEX08_06": 1.3908456310310637})
factors.update({"AVG_NDS_VOX_HEX08_12": 1.2743639355748966})
factors.update({"AVG_NDS_VOX_HEX08_24": 1.2596287755834055})
factors.update({"AVG_NDS_VOX_HEX08_48": 1.2179251259194124})
factors.update({"AVG_NDS_VOX_HEX08_72": 1.1793611129414954})
factors.update({"AVG_NDS_VOX_HEX08_96": 0.5021699727716449}) # not plausible!
#
# Solid FEs:
DOF = 3
#1025
DEBUG = False
#
#
USER_DEF_RATIO = 1.0
#
# List of optimal job sizes
# job_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
job_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
#
# Has to be open to lower counts as the topology (job_sizes) better has to be respected 
IDEAL_CCN = 32

# -----------------------------------------------------------------------------
# How many compute nodes are assigned to the domain. 1/64cn = 2c is already 
# accounted for.
# -----------------------------------------------------------------------------
core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072 ]
core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 ] # High ELEM COUNT
# core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192] 
# core_categories = [ 1024, 2048, 4096] # OOM prevention for size >= 72
# core_categories = [ 4096] # OOM prevention for size >= 96

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
parser.add_argument('-metaF', action="store", dest="meta_file", help="Input meta file.")
# parser.add_argument('-ncno', action="store", dest="ncno", help="Number of compute nodes.")
parser.add_argument('-size', action="store", dest="size", help="Size of the domains in format of, e.g., '12'.")

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

# if bool(cmd_args.ncno):
#     requested_job_size = int(cmd_args.ncno)
# else:
#     requested_job_size = 2

if bool(cmd_args.size):
    SIZE = str(cmd_args.size)
else:
    print("You need to define a domain size.")
    exit(3)

# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")
# size_dmn =  keywords.get("SIZE_DOMAIN")

basename = os.path.splitext(cmd_args.meta_file)[0]

vox_file    = basename + ".vox"
cores_file  = basename + ".cores"
meta_file     = basename + ".meta"

print("-- Outputfiles: ")
print("-- meta_file:   ", meta_file)
print("-- vox_file:    ", vox_file)
print("-- cores_file:  ", cores_file)
print("--")

# -----------------------------------------------------------------------------
# Read the information per domain
# -----------------------------------------------------------------------------
file_missing = False
if not os.path.exists(vox_file):
    print("EE Input file »" + vox_file + "« does not exist.")
    file_missing = True

if file_missing: 
    print("-- Program stopped")
    print(FMT_SEP)
    exit(2)

vox_list = ik8_to_list(vox_file)

print("-- Entries in vox_list:     ", len(vox_list))
print(FMT_SEP)

# -----------------------------------------------------------------------------
# Define suggested numbers of els per core
# -----------------------------------------------------------------------------
# suggested = 10000
from_baseline_analyses = [6659, 5814, 4148, 5622, 7419, 8487]
sizes_sug              = ["96", "72", "48", "24", "12", "06"]

suggested = 0
for ii in range(len(sizes_sug)):
    if SIZE == sizes_sug[ii]:
        suggested = from_baseline_analyses[ii]

if suggested < 2:
    print("EE Error while defining best no of elements per core.")
    print(FMT_SEP)
    exit(4)






print("-- Suggested voxel/FEs per core:", suggested) # Elements per part
print(FMT_SEP)

list_of_NoofComms = [ int(1) for ii in range(len(core_categories))]

# number of cores per domain
list_of_comms = core_categories

# number of domains of each communicator
list_of_NoDmns = [ int(0) for ii in range(len(core_categories))]

cores_list = []
for vox in vox_list:
    ideal_core_no = vox/suggested

    # lo_no_ideal_core = vox/suggested_lo
    # up_no_ideal_core = vox/suggested_up

    if ideal_core_no < min(core_categories):
        ideal_core_no = min(core_categories)
        
    if ideal_core_no < 2.0:
        ideal_core_no = 2

    for ii in range(len(core_categories)-1):

        if core_categories[ii] <= ideal_core_no and core_categories[ii+1] > ideal_core_no:

            cores_list.append(core_categories[ii])

            list_of_NoDmns[ii] += int(1)
            
            break

        if ii == (len(core_categories)-2):
            cores_list.append(core_categories[ii])
            elems = vox/core_categories[ii]
            print("WW Number of FEs per dmn :", elems)

# -----------------------------------------------------------------------------
# Refactor the list of the number of domains and of communicators
# -----------------------------------------------------------------------------
for ii in range(len(list_of_NoDmns)):
    if list_of_NoDmns[ii] == 0:
        list_of_NoofComms[ii] = 0

sum_of_cores = 0
for ii in range(len(list_of_comms)):
    sum_of_cores += list_of_comms[ii] * list_of_NoofComms[ii]

job_size = max(cores_list) / IDEAL_CCN

# -----------------------------------------------------------------------------
# User feedback
# -----------------------------------------------------------------------------
with open(meta_file, 'a') as f:
    f.write("w CORES_REQUESTED       " + str(sum_of_cores) + '\n')
f.close()

print("-- list_of_NoDmns:                      ", list_of_NoDmns)
print("-- list_of_comms:                       ", list_of_comms)
print("-- Sum of cores:                        ", sum_of_cores)
print("-- Number of compute nodes:             ", int(np.ceil(sum_of_cores/IDEAL_CCN)))


# -----------------------------------------------------------------------------
# Write the parts/cores per domain for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
f = open(cores_file, 'wb')
for ii in cores_list:
    f.write(int(ii).to_bytes(8, byteorder='little', signed=True))
f.close()
#

# -----------------------------------------------------------------------------
# Prepare provenance basename jobs
# Create corresponding meta files and directories
# -----------------------------------------------------------------------------
for ii in range(len(list_of_comms)):
    if list_of_comms[ii] < 0.1:
        continue

    if list_of_NoDmns[ii] == 0:
        continue

    new_basename = basename + "-" + str(list_of_comms[ii]) + "c"
    new_meta_file = new_basename + ".meta"
    shutil.copy(cmd_args.meta_file, new_meta_file)

    # -----------------------------------------------------------------------------
    # Set keywords of provenance files
    # -----------------------------------------------------------------------------  
    file_path = new_meta_file

    keyword_to_replace = 'MESH_PER_SUB_DMN'
    new_line_content = "r MESH_PER_SUB_DMN      " + str(list_of_comms[ii])
    replace_line_in_file(file_path, keyword_to_replace, new_line_content)

    keyword_to_replace = 'COMPUTE_NODES'
    new_line_content = "r COMPUTE_NODES         " + str(list_of_comms[ii]/IDEAL_CCN)
    replace_line_in_file(file_path, keyword_to_replace, new_line_content)


    # -----------------------------------------------------------------------------
    # Create directories now
    # -----------------------------------------------------------------------------
    main_process = 0
    worker_main_rank = 1

    path_spec = basename + "/" + f"Rank_{worker_main_rank:07}"
    Path(path_spec).mkdir(parents=True, exist_ok=True)

    no_cores_total = (job_size*IDEAL_CCN)+1

    remainder = no_cores_total % list_of_comms[ii]

    no_folders = int(no_cores_total - remainder)

    worker_main_rank = 1
    for jj in range(no_folders):

        path_spec = new_basename + "/" + f"Rank_{worker_main_rank:07}"
        Path(path_spec).mkdir(parents=True, exist_ok=True)

        worker_main_rank += list_of_comms[ii]


end_time = time.time()
#
elapsed_time = end_time - start_time
#


if len(vox_list) != len(cores_list):
    print()
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    print("EE *.cores and *.vox are not of equal length! DTC will stop. EE")
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    exit(10)
else:
    print(FMT_SEP)