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
# Cutting stocks --> Time "logs" of best length (shortest necessary) for 
# highest utilization of the cluster. 
#------------------------------------------------------------------------------
# DTC distributes the domains in a round-robin manner. Therefore, the idling
# time of the processors must be minimized, because the efficiency gain by 
# adjusting the comm size to the ideal cpds have to overcompensate the idling
# processors.
#------------------------------------------------------------------------------
import os, struct, argparse, sys, time
import numpy as np
import pandas as pd
from pathlib import Path
import struct
#
sys.path.insert(1, '/home/geb/00_bone_eval_chain/P_MOD_Python')
#
from mod_decomposition import *
from mod_mechanical import *
from mod_meta_parser import *
#
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
job_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
#
# Has to be open to lower counts as the topology (job_sizes) better has to be respected 
IDEAL_CCN = 32
# IDEAL_CCN = 16

# -----------------------------------------------------------------------------
# How many compute nodes are assigned to the domain. 1/64cn = 2c is already 
# accounted for.
# -----------------------------------------------------------------------------
core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072 ]
core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 ] # High ELEM COUNT
# core_categories = [ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
# core_categories = [  1024, 2048, 4096] # OOM prevention for size >= 72
# core_categories = [   2048, 4096] # OOM prevention for size >= 72

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
parser.add_argument('-ncno', action="store", dest="ncno", help="Number of compute nodes.")
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

if bool(cmd_args.ncno):
    requested_job_size = int(cmd_args.ncno)
else:
    requested_job_size = 2

if bool(cmd_args.size):
    SIZE = str(cmd_args.size)
else:
    print("You need to define a domain size.")
    exit(3)

# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")

mi_el_type  =  keywords.get("MICRO_ELMNT_TYPE")
ma_el_order =  keywords.get("MACRO_ELMNT_ORDER")
size_dmn    =  keywords.get("SIZE_DOMAIN")

basename = os.path.splitext(cmd_args.meta_file)[0]

vox_file    = basename + ".vox"
groups_file = basename + ".groups"
comms_file  = basename + ".comms"
cores_file  = basename + ".cores"
meta_file     = basename + ".meta"

print("-- Outputfiles: ")
print("-- meta_file:   ", meta_file)
print("-- vox_file:    ", vox_file)
print("-- comms_file:  ", comms_file)
print("-- cores_file:  ", cores_file)
print("-- groups_file: ", groups_file)
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
# 10000 is ok, 5000 still pretty good. This allows for a naive approach in 
# which the number of vox/core varies between 5000 and 10000 by rounding to
# max twice as many cores as required
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

    # if lo_no_ideal_core < 2.0:
    #     lo_no_ideal_core = 2
    # if up_no_ideal_core < 2.0:
    #     up_no_ideal_core = 2
    found_assignment = False
    for ii in range(len(core_categories)-1):

        if core_categories[ii] <= ideal_core_no and core_categories[ii+1] > ideal_core_no:

            cores_list.append(core_categories[ii])

            list_of_NoDmns[ii] += int(1)
            
            found_assignment = True

            break
        if ii == (len(core_categories)-2):
            cores_list.append(core_categories[ii])
            elems = vox/core_categories[ii]
            print("WW Number of FEs per dmn :", elems)

# -----------------------------------------------------------------------------
# Refactor the list of the number of communicators
# -----------------------------------------------------------------------------
for ii in range(len(list_of_NoDmns)):
    if list_of_NoDmns[ii] == 0:
        list_of_NoofComms[ii] = 0

# -----------------------------------------------------------------------------
# Refactor the list of the number of communicators
# -----------------------------------------------------------------------------
sum_of_cores = 0
for ii in range(len(list_of_comms)):
    sum_of_cores += list_of_comms[ii] * list_of_NoofComms[ii]

temp_job_size = int(np.ceil((sum_of_cores)/IDEAL_CCN))


if sum_of_cores < IDEAL_CCN:
    print("-- Initial sum_of_cores lower than IDEAL_CCN")
    print(FMT_SEP)
    optimize = True
else: 
    optimize = False


no_of_serial_domains = []

if requested_job_size == 0 and not optimize:
    pass
elif temp_job_size > requested_job_size:
    print("WW Requested job_size smaller than optimal one.")
    print(FMT_SEP)
else:

    go_on = True
    while go_on:              
        

        # Get list of domains computed serially (which relates to the walltime)
        for ii in range(len(list_of_NoDmns)):

            if list_of_NoofComms[ii] == 0:
                no_of_serial_domains.append(0)
            else:    
                no_of_serial_domains.append(list_of_NoDmns[ii]/list_of_NoofComms[ii])

        max_serial = max(no_of_serial_domains)

        # Reduce the worst one
        for ii in range(len(no_of_serial_domains)):
            if max_serial == no_of_serial_domains[ii]:
                list_of_NoofComms[ii] += int(1)

        # Get the new sum of cores required
        sum_of_cores = 0
        for ii in range(len(list_of_comms)):
            sum_of_cores += list_of_comms[ii] * list_of_NoofComms[ii]

        # Check the job size and stop if optimal
        temp_job_size = int(np.ceil((sum_of_cores)/IDEAL_CCN))

        if temp_job_size > requested_job_size or temp_job_size == max(job_sizes):
            print("WW Requested job size not met.")
            print(FMT_SEP)
            go_on = False
        elif temp_job_size == requested_job_size:
            go_on = False


job_size = temp_job_size 
sum_of_cores += 1 # To account for the main rank

# -----------------------------------------------------------------------------
# Write the parts/cores per domain for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
no_of_domains_checksum = 0
last_value = 0
f = open(groups_file, 'wb')
for ii in range(len(list_of_comms)):
    if ii > 0 and last_value == list_of_comms[ii]:
        continue

    f.write(list_of_comms[ii].to_bytes(8, byteorder='little', signed=True))
    f.write(list_of_NoofComms[ii].to_bytes(8, byteorder='little', signed=True))
    f.write(list_of_NoDmns[ii].to_bytes(8, byteorder='little', signed=True))
    
    no_of_domains_checksum += list_of_NoDmns[ii]

    last_value = list_of_comms[ii]

f.close()



# -----------------------------------------------------------------------------
# User feedback
# -----------------------------------------------------------------------------
with open(meta_file, 'a') as f:
    f.write("w CORES_REQUESTED       " + str(sum_of_cores) + '\n')
    f.write("w NCNO                  " + str(job_size) + '\n')
f.close()


print("-- Number of domains written to files:  ", no_of_domains_checksum)
print("-- list_of_NoDmns:                      ", list_of_NoDmns)
print("-- list_of_comms:                       ", list_of_comms)
print("-- list_of_NoofComms:                   ", list_of_NoofComms)
if no_of_serial_domains:
    print("-- Number of serial domains:            ", np.around(no_of_serial_domains))
print("-- Sum of cores:                        ", sum_of_cores)
print("-- Number of compute nodes:             ", job_size)

# -----------------------------------------------------------------------------
# Get an impression on the wasted core hours
# -----------------------------------------------------------------------------
# if no_of_serial_domains:
#     max_serial = max(no_of_serial_domains)

#     tot_core_h_std_units = 0
#     wasted_instances = 0

#     print(no_of_serial_domains)
#     print(list_of_NoofComms)
#     print(list_of_comms)

#     print(list_of_NoDmns)
#     print(list_of_comms)

#     print(FMT_SEP)

#     for ii in range(len(no_of_serial_domains)):
#         tot_core_h_std_units += (list_of_NoDmns[ii] * list_of_comms[ii])
#         wasted_instances += ((max_serial - no_of_serial_domains[ii]) * list_of_NoofComms[ii] * list_of_comms[ii])

#     print(FMT_SEP)
#     print("-- Max number of serial domains:        ", int(max_serial))
#     print("-- Total core hours (std time units):   ", tot_core_h_std_units)
#     print("-- 'Wasted' core hors (std time units): ", int(round(wasted_instances,0)))
#     print("-- Percentage 'wasted':                 ", round(wasted_instances/tot_core_h_std_units*100,3), "%")
#     print(FMT_SEP)

# -----------------------------------------------------------------------------
# Write the parts/cores per domain for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
f = open(cores_file, 'wb')
for ii in cores_list:
    f.write(int(ii).to_bytes(8, byteorder='little', signed=True))
f.close()
#

# -----------------------------------------------------------------------------
# Write the comms file for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
f = open(comms_file, 'wb')
for ii in range(len(cores_list)):
    f.write((0).to_bytes(4, byteorder='little', signed=True))
f.close()


# -----------------------------------------------------------------------------
# Create directories now
# -----------------------------------------------------------------------------
main_process = 0
worker_main_rank = 1

path_spec = basename + "/" + f"Rank_{worker_main_rank:07}"
Path(path_spec).mkdir(parents=True, exist_ok=True)

for ii in range(len(list_of_comms)):

    cpd = list_of_comms[ii]
    comms_of_cpd = list_of_NoofComms[ii]

    for jj in range(comms_of_cpd):
        worker_main_rank += cpd

        path_spec = basename + "/" + f"Rank_{worker_main_rank:07}"
        Path(path_spec).mkdir(parents=True, exist_ok=True)

# Simple workaround to remove the last directory, which is not needed.
Path(path_spec).rmdir()

path_spec = basename + "/results_domains"
Path(path_spec).mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Create directories now
# -----------------------------------------------------------------------------
end_time = time.time()
#
elapsed_time = end_time - start_time
#
print("-- Single core runtime of this script:  ", round(elapsed_time,1), "s")

if len(vox_list) != len(cores_list):
    print()
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    print("EE *.cores and *.vox are not of equal length! DTC will stop. EE")
    print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    exit(10)
else:
    print(FMT_SEP)