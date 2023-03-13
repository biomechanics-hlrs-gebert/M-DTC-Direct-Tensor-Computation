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
import pandas as pd
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
# Solid FEs:
DOF = 3
#
# List of optimal job sizes
job_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
#
# Target parts per compute node [min, optimal, max]
# target_pcn = [15, 20, 25]
#
# Has to be open to lower counts as the topology (job_sizes) better has to be respected 
target_pcn = [10, 20, 22]
#
range_of_FEs_per_part = [1500,4000,500]
FEs_part = range_of_FEs_per_part
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
# ChatGPT:
def cutting_stock(item_lengths, stock_length):
    # sort the items in decreasing order
    item_lengths.sort(reverse=True)
    
    # initialize the number of bins used and the current bin capacity
    num_bins = 0
    current_bin_capacity = 0
    
    # iterate through the items
    for item_length in item_lengths:
        # try to fit the item in the current bin
        if current_bin_capacity >= item_length:
            current_bin_capacity -= item_length
        else:
            # if the item doesn't fit in the current bin, start a new bin
            num_bins += 1
            current_bin_capacity = stock_length - item_length
    
    return num_bins
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
parser.add_argument('-ppn', action="store", dest="ppn", help="Input meta file).")
parser.add_argument('-ncno', action="store", dest="ncno", help="Number of compute nodes).")

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
#
no_of_dmns_in_comm = 0
runtime_in_comm = 0
#
# Derived with analyze_runtime_to_BVTV
dmn_size_avg=sum([float(i) for i in size_dmn])/3
#
if dmn_size_avg == 0.6:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_06")
elif dmn_size_avg == 1.2:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_12")
elif dmn_size_avg == 2.4:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_24")
elif dmn_size_avg == 4.8:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_48")
elif dmn_size_avg == 7.2:
    node_factor   = factors.get("AVG_NDS_VOX_HEX08_72")

# -----------------------------------------------------------------------------
# This is a normalized measure of the expected runtime. All domains use 
# the same macro element order, the same micro element type and about the same
# number of finite element nodes per domain. Consequently, the runtime 
# expected is roughly the same for all. The problem collapses to a much simpler
# version of the cutting stock problem.
#
# The problem posed will become a »true« cutting stock problem if the runtimes
# change the number of FE nodes per domain. But this max only happen in a 
# latter step.
# -----------------------------------------------------------------------------
# The runtime however is not removed from this code to allow for such 
# optimization.
# -----------------------------------------------------------------------------
total_factors = tf_nds * tf_ma_el 

#
# -----------------------------------------------------------------------------
# Loop over FE nodes per domain for optimizing the topology aware batch 
# scheduling
# -----------------------------------------------------------------------------
for FE_nodes_part in range(FEs_part[0],FEs_part[1]+1,FEs_part[2]):

    # Reset dataframes
    catalogued_data = pd.DataFrame(columns=['domain', 'ppd', 'expected runtime'])
    bins_per_ppd = pd.DataFrame(columns=['ppd', 'bins'])

    for dmn in range(len(dmn_no_list)):
        FE_nodes_dmn = vox_list[dmn] * node_factor

        parts_per_domain = int(np.floor(FE_nodes_dmn/FE_nodes_part))
        
        # expected runtime may be calibrated in a latter step
        expected_runtime = 1.0 * total_factors

        new_entry = pd.DataFrame({'domain': dmn, 'ppd': parts_per_domain, 'expected runtime': expected_runtime}, index=[0])
        catalogued_data = pd.concat([catalogued_data, new_entry])


    min_total_time_units = 9999999999999
    max_total_time_units = 0
    unique_ppds = catalogued_data['ppd'].unique()

    # Batches of unique parts per domain
    for uppd in unique_ppds:

        df_ppd = catalogued_data[catalogued_data['ppd'] == uppd]

        expRuntime = df_ppd['expected runtime'].sum()

        # print("-- No of domains with " + str(uppd) + " parts:", len(df_ppd), \
        #     "with a total runtime of", expRuntime, "std time units")

        if min_total_time_units > expRuntime:
            min_total_time_units = expRuntime
        if max_total_time_units < expRuntime:
            max_total_time_units = expRuntime

    target_time_budget = min_total_time_units

    for uppd in unique_ppds:

        df_dmns = catalogued_data[catalogued_data['ppd'] == uppd]

        num_bins = cutting_stock(df_dmns["expected runtime"].tolist(), target_time_budget)

        new_entry = pd.DataFrame({'ppd': uppd, 'bins': num_bins}, index=[0])
        bins_per_ppd = pd.concat([bins_per_ppd, new_entry])

    # Number of cores required with this setup
    sum_of_cores = sum(bins_per_ppd["ppd"]*bins_per_ppd["bins"])

    # -----------------------------------------------------------------------------
    # Search for optimum packaging with topology aware compute node numbers
    # -----------------------------------------------------------------------------
    min_delta_cores = 999999999
    curr_job_size = 0
    curr_target_pcn = 0
    for js in job_sizes:
        for ii in range(target_pcn[0], target_pcn[2]):
            cores_avail = ii*js

            delta_cores = cores_avail-sum_of_cores

            if delta_cores < min_delta_cores and delta_cores >= 0:
                min_delta_cores = delta_cores

                curr_job_size = js
                curr_target_pcn = ii

    sum_of_cores = sum(bins_per_ppd["ppd"]*bins_per_ppd["bins"])
    print("--", bins_per_ppd)
    print("-- FE_nodes_part:   ", FE_nodes_part)
    print("-- sum_of_cores:   ", sum_of_cores)
    print("-- curr_job_size:  ", curr_job_size)
    print("-- curr_target_pcn:", curr_target_pcn)
    print("-- ")
    print("-- Cores available:", curr_job_size*curr_target_pcn)
    print(FMT_STRING)


# # -----------------------------------------------------------------------------
# # Retrieve optimal solution
# # -----------------------------------------------------------------------------
# catalogued_data = pd.DataFrame(columns=['domain', 'ppd', 'expected runtime'])
# bins_per_ppd = pd.DataFrame(columns=['ppd', 'bins'])
# for dmn in range(len(dmn_no_list)):
#     FE_nodes_dmn = vox_list[dmn] * node_factor

#     parts_per_domain = int(np.floor(FE_nodes_dmn/best_FE_nodes_part))
    
#     # expected runtime may be calibrated in a latter step
#     expected_runtime = 1.0 * total_factors

#     new_entry = pd.DataFrame({'domain': dmn, 'ppd': parts_per_domain, 'expected runtime': expected_runtime}, index=[0])
#     catalogued_data = pd.concat([catalogued_data, new_entry])


# min_total_time_units = 9999999999999
# max_total_time_units = 0
# unique_ppds = catalogued_data['ppd'].unique()

# # Batches of unique parts per domain
# for uppd in unique_ppds:

#     df_ppd = catalogued_data[catalogued_data['ppd'] == uppd]

#     expRuntime = df_ppd['expected runtime'].sum()

#     print("-- No of domains with " + str(uppd) + " parts:", len(df_ppd), \
#         "with a total runtime of", expRuntime, "std time units")

#     if min_total_time_units > expRuntime:
#         min_total_time_units = expRuntime
#     if max_total_time_units < expRuntime:
#         max_total_time_units = expRuntime

# target_time_budget = min_total_time_units

# for uppd in unique_ppds:

#     df_dmns = catalogued_data[catalogued_data['ppd'] == uppd]

#     num_bins = cutting_stock(df_dmns["expected runtime"].tolist(), target_time_budget)

#     new_entry = pd.DataFrame({'ppd': uppd, 'bins': num_bins}, index=[0])
#     bins_per_ppd = pd.concat([bins_per_ppd, new_entry])

# # Number of cores required with this setup
# sum_of_cores = sum(bins_per_ppd["ppd"]*bins_per_ppd["bins"])

# print(FMT_STRING)
# print("--", bins_per_ppd)
# print(FMT_STRING)

# print("-- sum_of_cores:   ", sum_of_cores)
# print("-- best_job_size:  ", best_job_size)
# print("-- best_target_pcn:", best_target_pcn)
# print("-- ")
# print("-- Cores available:", best_job_size*best_target_pcn)

