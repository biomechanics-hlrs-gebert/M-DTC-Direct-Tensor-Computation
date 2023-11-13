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
# Has to be open to lower counts as the topology (job_sizes) better has to be respected 
IDEAL_CCN = 32


# -----------------------------------------------------------------------------
# How many compute nodes are assigned to the domain. 1/64cn = 2c is already 
# accounted for.
# -----------------------------------------------------------------------------
core_cat = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072 ]
core_cat = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 ] # High ELEM COUNT
# core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192] 
# core_categories = [ 1024, 2048, 4096] # OOM prevention for size >= 72
# core_categories = [ 4096, 8192 ] # OOM prevention for size >= 96


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
parser.add_argument('-cn', action="store", dest="cns", help="Number of compute nodes.")
parser.add_argument('-size', action="store", dest="size", help="Set of size 06, 12, 24, 48.")

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

if not bool(cmd_args.cns):
    print("Please specify the number of compute nodes.")
    exit(1)
else:
    CN = int(cmd_args.cns)    

if not bool(cmd_args.size):
    print("Please specify the size of the set (06, 12, 24, 48).")
    exit(1)
else:
    size = cmd_args.size

if not os.path.exists(cmd_args.meta_file):
    print("Input file »" + cmd_args.meta_file + "« does not exist.")
    exit(2)



# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")

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
core_categories = []
for ii in range(len(core_cat)):
    if core_cat[ii] > IDEAL_CCN*CN:
        break
    core_categories.append(core_cat[ii])


nel_core_min =  500
nel_core_max = 35000

if size == "06":
    TGT_CORES = 16
elif size == "12":
    TGT_CORES = 32
elif size == "24":
    TGT_CORES = 64
elif size == "48":
    TGT_CORES = 256

for ii in range(len(core_categories)):
    if core_categories[ii] == TGT_CORES:
        TGT_II = ii

list_of_NoofComms = [ int(1) for ii in range(len(core_categories))]

# number of cores per domain
list_of_comms = core_categories

# number of domains of each communicator
list_of_NoDmns = [ int(0) for ii in range(len(core_categories))]

par_dmns = []
cores_list = []
for vox in vox_list:
    tc = -1
    el_number = int(np.around(vox/TGT_CORES,0))

    if el_number < nel_core_min and el_number > 1:

        # for tgt_cores in reversed([1,2,4,8,16]):
        for tgt_cores in reversed([ 2**ii  for ii in range(int(np.log2(TGT_CORES)))]):
            
            el_number = int(np.around(vox/tgt_cores,0))

            if el_number > nel_core_min:
                tc = tgt_cores
                break

    elif el_number > nel_core_max:

        # for tgt_cores in reversed([64,128,256,512,1024,2048,4096]):
        for tgt_cores in reversed([ 2**ii for ii in range(int(np.log2(TGT_CORES*2)) , int(np.log2(CN*IDEAL_CCN*2)))]):
            
            el_number = int(np.around(vox/tgt_cores,0))

            if el_number < nel_core_max:
                tc = tgt_cores
                break
        # Try the largest number of cores available
        if tc == -1:
            tc = CN*IDEAL_CCN

    elif el_number < 1:

        tc = 0

    else:
        list_of_NoDmns[TGT_II] += 1
        tc = TGT_CORES


    if tc < 2:
        tc = 2

    cores_list.append(tc)

    # -----------------------------------------------------------------------------
    # Refactor the list of the number of domains and of communicators
    # -----------------------------------------------------------------------------
    if core_categories[TGT_II] != tc:

        for ii in range(len(core_categories)):

            if core_categories[ii] == tc:
                list_of_NoDmns[ii] += 1
                break

    # print("WW Number of FEs per dmn :", elems)

# -----------------------------------------------------------------------------
# If no of domains of a comm size below 2 --> Do not even start the computation
# Simply delete the domains (not too many lost for a perf investigation)
# -----------------------------------------------------------------------------
list_of_NoDmns = [ 0 if ii<2 else ii  for ii in list_of_NoDmns ] 

# -----------------------------------------------------------------------------
# Reduce number of cores if not enough domains requested. Leads to idling CNs
# -----------------------------------------------------------------------------
par_dmns = [ (list_of_NoDmns[ii]*list_of_NoofComms[ii])/(CN*IDEAL_CCN) for ii in range(len(list_of_NoDmns)) ]
cores_requested = []
for ii in range(len(par_dmns)):
    if par_dmns[ii] < 2.0:

        if list_of_NoDmns[ii] == 0:
            cores_requested.append(0)
        else:
            c_req = int(np.ceil(list_of_NoDmns[ii]/2) * list_of_comms[ii])+1

            if c_req > CN*TGT_CORES+1:
                cores_requested.append(CN*IDEAL_CCN+1)
            else:
                cores_requested.append(c_req)

            par_dmns[ii] = int((cores_requested[ii]-1) / list_of_comms[ii])
    else:
        cores_requested.append(CN*TGT_CORES+1)

# -----------------------------------------------------------------------------
# Refactor the list of the number of domains and of communicators
# -----------------------------------------------------------------------------
for ii in range(len(list_of_NoDmns)):
    if list_of_NoDmns[ii] == 0:
        list_of_NoofComms[ii] = 0

job_size = max(cores_requested)

print("-- core_categories:                     ", core_categories)
print("-- list_of_NoDmns:                      ", list_of_NoDmns)
print("-- list_of_comms:                       ", list_of_comms)
print("-- parallel domains:                    ", np.around(par_dmns,2))
print("-- Cores requested:                     ", cores_requested)
print("-- Total domains scheduled:             ", sum(list_of_NoDmns))
print("-- Target number of cores per domain:   ", TGT_CORES)
print("-- Job size:                            ", job_size)
print("-- No of cn, checksum:                  ", (job_size-1)/IDEAL_CCN)
print("-- Min number of elements per core:     ", nel_core_min)
print("-- Max number of elements per core:     ", nel_core_max)

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

    remainder = cores_requested[ii] % list_of_comms[ii]

    no_folders = int((cores_requested[ii] - remainder)/list_of_comms[ii])

    worker_main_rank = 1
    for jj in range(no_folders):

        path_spec = new_basename + "/" + f"Rank_{worker_main_rank:07}"
        Path(path_spec).mkdir(parents=True, exist_ok=True)

        worker_main_rank += list_of_comms[ii]


end_time = time.time()
#
elapsed_time = end_time - start_time

if len(vox_list) != len(cores_list):
    print()
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    print("EE *.cores and *.vox are not of equal length! DTC will stop. EE")
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    exit(10)
else:
    print(FMT_SEP)

print()