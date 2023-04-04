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
# adjusting the comm size to the ideal ppds have to overcompensate the idling
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
#
DEBUG = False
#
#
USER_DEF_RATIO = 1.0
#
# List of optimal job sizes
job_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
#
# Target parts per compute node [min, optimal, max]
# target_pcn = [15, 20, 25]
#
# Has to be open to lower counts as the topology (job_sizes) better has to be respected 
target_pcn = [5, 20, 100]
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
def hist(vox, min, max, no_bins):

    bins = np.linspace(min, max, no_bins)
    weightsa = np.ones_like(vox)/float(len(vox))
    histogram_data = np.histogram(np.array(vox), bins, weights = weightsa)

    return(histogram_data)
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
FMT_SEP="-------------------------------------------------------------------------------"
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
# parser.add_argument('-ppn', action="store", dest="ppn", help="Input meta file).")
parser.add_argument('-ncno', action="store", dest="ncno", help="Number of compute nodes.")
parser.add_argument('-bins', action="store", dest="bins", help="Number of different sized of MPI groups.")
parser.add_argument('-idle-offset', action="store", dest="idle_offset", help="Offset of min runtime to target runtime (0.1 => 1.1*min runtime = target runtime)")
parser.add_argument('-debug', action="store", dest="debug", help="Number of compute nodes.")

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

# if bool(cmd_args.ppn):
#     print("MM Optimizing with user defined number of parts per compute node.")
#     target_pcn = [int(cmd_args.ppn)]

if bool(cmd_args.ncno):
    job_sizes = [int(cmd_args.ncno)]


if bool(cmd_args.debug):
    if cmd_args.debug[0].lower() == "y":
        DEBUG = True

    print("MM User switched on DEBUG = True.")

if not bool(cmd_args.idle_offset):
    idle_offset = 0.001
else:
    idle_offset = float(cmd_args.idle_offset)

if not bool(cmd_args.bins):
    bins = list(range(3, 16))
else:
    bins = [int(cmd_args.bins)]


# -----------------------------------------------------------------------------
# Parse the meta file
# -----------------------------------------------------------------------------
keywords = parse_meta(cmd_args.meta_file, "TENSOR_COMPUTATION")

mi_el_type  =  keywords.get("MICRO_ELMNT_TYPE")
ma_el_order =  keywords.get("MACRO_ELMNT_ORDER")
size_dmn    =  keywords.get("SIZE_DOMAIN")

basename = os.path.splitext(cmd_args.meta_file)[0]

dmn_no_file = basename + ".status"
vox_file    = basename + ".vox"
groups_file = basename + ".groups"
comms_file  = basename + ".comms"
parts_file  = basename + ".parts"
exp_rtex_file = basename + ".rtex"
exp_rtac_file = basename + ".rtac"
meta_file     = basename + ".meta"

print("-- Outputfiles: ")
print("-- meta_file:   ", meta_file)
print("-- vox_file:    ", vox_file)
print("-- comms_file:  ", comms_file)
print("-- parts_file:  ", parts_file)
print("-- dmn_no_file: ", dmn_no_file)
print("-- groups_file: ", groups_file)
print("-- exp_rtex_file: ", exp_rtex_file)
print("-- exp_rtac_file: ", exp_rtac_file)
print("--")

# -----------------------------------------------------------------------------
# Read the information per domain
# -----------------------------------------------------------------------------
file_missing = False
if not os.path.exists(vox_file):
    print("EE Input file »" + vox_file + "« does not exist.")
    file_missing = True

if not os.path.exists(dmn_no_file):
    print("EE Input file »" + dmn_no_file + "« does not exist.")
    file_missing = True

if file_missing: 
    print("-- Program stopped")
    print(FMT_SEP)
    exit(2)

nds_list = ik8_to_list(vox_file)
dmn_no_list = ik8_to_list(dmn_no_file)

print("-- Entries in nds_list:     ", len(nds_list))
print("-- Entries in dmn_no_list:  ", len(dmn_no_list))
print(FMT_SEP)

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
print(FMT_SEP)

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
    node_factor = factors.get("AVG_NDS_VOX_HEX08_06")
    suggested = 3000
elif dmn_size_avg == 1.2:
    node_factor = factors.get("AVG_NDS_VOX_HEX08_12")
    suggested = 6000
elif dmn_size_avg == 2.4:
    node_factor = factors.get("AVG_NDS_VOX_HEX08_24")
    suggested = 6000
elif dmn_size_avg == 4.8:
    node_factor = factors.get("AVG_NDS_VOX_HEX08_48")
    suggested = 6000
elif dmn_size_avg == 7.2:
    node_factor = factors.get("AVG_NDS_VOX_HEX08_72")
    suggested = 6000
elif dmn_size_avg == 9.6:
    node_factor = factors.get("AVG_NDS_VOX_HEX08_72")
    suggested = 6000

# Currently not in use
FEs_part[0] = suggested / 2 
FEs_part[1] = suggested * 2 
FEs_part[2] = suggested / 8

print("-- Suggested pcno:", suggested)
print(FMT_SEP)

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
best_wct = 999999999
best_job_size = 0
best_target_pcn = 0
best_sum_of_cores = 0

# print("min(FE_nodes_dmn)", min(FE_nodes_dmn))
# print("max(FE_nodes_dmn)", max(FE_nodes_dmn))
# print("-- sum(b_a[0]):" ,sum(bins_available[0]))
# print("-- sum(b_a[1]):" ,sum(bins_available[1]))
# print("-- bins:" ,bins_available)
# exit(0)

# for FE_nodes_part in range(FEs_part[1],FEs_part[1]+1,FEs_part[2]):

# -----------------------------------------------------------------------------
# Adjust the list of voxels per domain by the node_factor. Otherwise, the 
# histogram will fail to cover all domains.
# -----------------------------------------------------------------------------
for ii in range(len(nds_list)):
    nds_list[ii] *= node_factor

for bin in bins:

    # Reset dataframes
    catalogued_data = pd.DataFrame(columns=['domain', 'ppd', 'expected runtime'])
    comms_per_ppd = pd.DataFrame(columns=['ppd', 'groups', 'No of Domains'])

    # Create the histogram data
    bins_avail = hist(nds_list, min(nds_list), max(nds_list), bin)
    bins_avail = bins_avail[1]

    # calculate the parts per domain of a corresponding histogram bin
    ppds_avail = []
    time_budget_allocated = []
    time_budget_required = []
    for ii in range(len(bins_avail)-1):
        avg_FEs_in_hist_bin = (bins_avail[ii] + bins_avail[ii+1]) / 2

        ppd = int(avg_FEs_in_hist_bin/suggested)

        # DTC only accepts at least two parts per domain
        if ppd == 1:
            ppd += 1

        # if ppd % 2 != 0:
        #     ppd += 1

        ppds_avail.append(ppd)
        time_budget_allocated.append(0)
        time_budget_required.append(0)

    checksum = 0
    checksummm = 0
    # Assign the ppd to the domain
    for dmn in range(len(dmn_no_list)):

        FE_nodes_dmn = nds_list[dmn]
    
        checksummm += 1
        parts_per_domain = 0
        for ii in range(len(bins_avail)-1):

            # -1 and +1 because >= and <= will not work.
            if FE_nodes_dmn > bins_avail[ii]-1 and FE_nodes_dmn < bins_avail[ii+1]+1:
                parts_per_domain = ppds_avail[ii]
                checksum += 1

        # -----------------------------------------------------------------------------
        # Get and take the effective number of FE nodes per part into account
        # -----------------------------------------------------------------------------


                            # ###############################################################################
                            # ##### Inject target factors here ##############################################
                            # tf    = [1.15356, 1.77793, 2.4023,  2.82787, 3.25344, 3.04901, 2.84457, 2.13608, 1.42758, 1.59531, 1.76303, 2.54851, 3.33399, 3.69061]
                            # Parts = [ 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17]
                            # ###############################################################################
                            # ###############################################################################
                            # if parts_per_domain in Parts:
                            #     for pp in range(len(Parts)): 
                            #         if parts_per_domain == Parts[pp]:
                            #             USER_DEF_RATIO = tf[pp]

        if parts_per_domain == 0:
            delta_percentage = 0.0
            delta_ratio = 0.0
        else:
            eff_nodes_part = int(FE_nodes_dmn / parts_per_domain)
            delta_ratio = (suggested - eff_nodes_part) / suggested * USER_DEF_RATIO 


        # expected runtime may be calibrated in a latter step
        expected_runtime = 1.0 * total_factors / (1+(delta_ratio))
        # expected_runtime = 1.0 * total_factors * (1+(delta_ratio/1000))
        # expected_runtime = 1.0 * total_factors * delta_ratio

        new_entry = pd.DataFrame({'domain': dmn, 'ppd': parts_per_domain, 'expected runtime': expected_runtime}, index=[0])
        catalogued_data = pd.concat([catalogued_data, new_entry])

    # Search for bin specific runtimes
    min_tot_t = 9999999999999
    max_tot_t = 0
    for ii in range(len(ppds_avail)):
        uupd = ppds_avail[ii]

        df_ppd = catalogued_data[catalogued_data['ppd'] == uupd]
        expRuntime = df_ppd['expected runtime'].sum()
        time_budget_required[ii] = expRuntime * uupd # core seconds

        if min_tot_t > expRuntime:
            min_tot_t = expRuntime
        if max_tot_t < expRuntime:
            max_tot_t = expRuntime


    print("-- min_tot_t", round(min_tot_t,1), "    max_tot_t", round(max_tot_t,1))

    # The target_time_budget is a global one. Therefore, we have to loop over 
    # ppds_avail twice (once to get the global target and then to operate on it)
    target_time_budget = min_tot_t + (max_tot_t-min_tot_t) * idle_offset

    # Calculate the actual time budgets of all groups via the cutting stock problem.
    for ii in range(len(ppds_avail)):
        uupd = ppds_avail[ii]

        df_dmns = catalogued_data[catalogued_data['ppd'] == uupd]

        groups = cutting_stock(df_dmns["expected runtime"].tolist(), target_time_budget)

        # core seconds = number of MPI groups * size of MPI groups * target time budget
        # +1 to acco
        time_budget_allocated[ii] = (groups * uupd + 1) * target_time_budget

        # Not skipping empty communicators results in processors that are idling all the time.
        if groups == 0:
            continue

        new_entry = pd.DataFrame({'ppd': uupd, 'groups': groups, 'No of Domains': df_dmns.shape[0]}, index=[0])
        comms_per_ppd = pd.concat([comms_per_ppd, new_entry])

    
    # -----------------------------------------------------------------------------
    # Number of cores required with this setup. +1 to account for the master 
    # process in the monolithic topology aware approach.
    # -----------------------------------------------------------------------------
    sum_of_cores = sum(comms_per_ppd["ppd"]*comms_per_ppd["groups"]) + 1

    # -----------------------------------------------------------------------------
    # Search for optimum packaging with topology aware compute node numbers
    # -----------------------------------------------------------------------------
    wct = 9999999999999999
    curr_wct = 9999999999999999
    curr_job_size = 0
    curr_target_pcn = 0
    ttb_allocated =  sum(time_budget_allocated) # total time budget (ttb)
    ttb_required =  sum(time_budget_required)

    for js in job_sizes:
        for ii in range(target_pcn[0], target_pcn[2]):
            cores_avail = ii*js

            delta_cores = cores_avail-sum_of_cores

            idle_time_delta_cores = delta_cores * target_time_budget

            # Allocated = all 
            total_idle_core_time = ttb_allocated - ttb_required + idle_time_delta_cores

            # delta cores negative --> more cores needed than available
            if total_idle_core_time < wct and delta_cores >= 0:
                wct = total_idle_core_time

                curr_job_size = js
                curr_target_pcn = ii
                curr_sum_of_cores = sum_of_cores
                curr_no_bins = bin
                curr_wct = wct # wasted core time
                curr_ttb_required = ttb_required
                curr_comms_per_ppd = comms_per_ppd
                curr_catalogued_data = catalogued_data
            else:
                continue

    if DEBUG:
        sum_of_cores = sum(curr_comms_per_ppd["ppd"]*curr_comms_per_ppd["groups"])
        print("--", curr_comms_per_ppd)
        print("-- Number of bins: ", curr_no_bins)
        print("-- sum_of_cores:   ", sum_of_cores)
        print("-- Cores available:", curr_job_size*curr_target_pcn)
        print("-- curr_job_size:  ", curr_job_size)
        print("-- curr_target_pcn:", curr_target_pcn)
        print("-- ")
        print("-- curr_wct:", curr_wct)
        print("-- curr_ttb_required:", curr_ttb_required)
        print(FMT_SEP)

    if curr_wct < best_wct and curr_wct >= 0:
        best_wct = curr_wct

        best_job_size = curr_job_size
        best_target_pcn = curr_target_pcn
        best_sum_of_cores = curr_sum_of_cores
        best_no_bins = curr_no_bins
        best_wct = wct
        best_ttb_required = ttb_required
        best_comms_per_ppd = curr_comms_per_ppd
        best_catalogued_data = curr_catalogued_data

print(FMT_SEP)

# -----------------------------------------------------------------------------
# User Feedback
# -----------------------------------------------------------------------------
fail = False
try:
    bins_used = hist(nds_list, min(nds_list), max(nds_list), best_no_bins)

    print("-- Cores available:        ", best_job_size*best_target_pcn)
    print("-- best_sum_of_cores:      ", best_sum_of_cores)
    print("-- Cores not in use:       ", (best_job_size*best_target_pcn)-best_sum_of_cores)
    print("-- best_job_size:          ", best_job_size, "  compute nodes")
    print("-- best_target_pcn:        ", best_target_pcn, "  parts per compute node")
    print("-- ")
    print("-- Number of bins:         ", best_no_bins)
    print("-- Bins used:              ", np.around(bins_used[1]))
    print("-- Core time idling/wasted:", round(best_wct,2 ), " Standard time units")
    print("-- Core time required:     ", round(best_ttb_required,2 ), " Standard time units")

    theo_eff = round(best_ttb_required/(best_ttb_required+best_wct)*100,3 )
    print("-- Theoretical efficiency: ", theo_eff, "%")
    print("--", best_comms_per_ppd)
except:
    fail = True

if best_sum_of_cores == 0 and bool(cmd_args.ncno) or fail:
    print("EE No solution found. User defined parameters may not fit.")
    print(FMT_SEP)
    exit(5)




if bool(cmd_args.ncno):
    print("MM Optimized with user defined number of compute nodes.")

print(FMT_SEP)



# -----------------------------------------------------------------------------
# Write MPI communicators to a binary file, which can be read easily by Fortran
# -----------------------------------------------------------------------------
list_of_comms = best_comms_per_ppd['ppd'].tolist()
list_of_NoofComms = best_comms_per_ppd['groups'].tolist()
list_of_NoDmns = best_comms_per_ppd['No of Domains'].tolist()

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

print("-- Number of domains written to files:", no_of_domains_checksum)

# -----------------------------------------------------------------------------
# Write the parts per domain for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
parts_list = best_catalogued_data['ppd'].tolist()
#
f = open(parts_file, 'wb')
for ii in parts_list:
    f.write((ii).to_bytes(8, byteorder='little', signed=True))
f.close()
#

# -----------------------------------------------------------------------------
# Write the expected runtimes to a binary file
# -----------------------------------------------------------------------------
exp_rt_list = best_catalogued_data['expected runtime'].tolist()
#
f = open(exp_rtex_file, 'wb')
for ii in exp_rt_list:
    float_num = struct.pack("f", ii)
    f.write(float_num)
f.close()

# -----------------------------------------------------------------------------
# Prepare the actual runtimes for use by DTC
# -----------------------------------------------------------------------------
exp_rt_list = best_catalogued_data['expected runtime'].tolist()
#
float_num = struct.pack("f", 0.0)
f = open(exp_rtac_file, 'wb')
for ii in range(len(parts_list)):
    f.write((0).to_bytes(8, byteorder='little', signed=True))
f.close()

# -----------------------------------------------------------------------------
# Write the comms file for reading by Fortran to a binary file
# -----------------------------------------------------------------------------
f = open(comms_file, 'wb')
for ii in range(len(parts_list)):
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

    ppd = list_of_comms[ii]
    comms_of_ppd = list_of_NoofComms[ii]

    for jj in range(comms_of_ppd):
        worker_main_rank += ppd

        path_spec = basename + "/" + f"Rank_{worker_main_rank:07}"
        Path(path_spec).mkdir(parents=True, exist_ok=True)

# Simple workaround to remove the last directory, which is not needed.
Path(path_spec).rmdir()

path_spec = basename + "/results_domains"
Path(path_spec).mkdir(parents=True, exist_ok=True)



# -----------------------------------------------------------------------------
# Update the meta file with the expected efficiency
# -----------------------------------------------------------------------------
with open(meta_file, 'a') as f:
    f.write("w EXPECTED_EFFICIENCY   " + str(theo_eff/100.0) + '\n')
    f.write("w CORES_REQUESTED       " + str(best_sum_of_cores) + '\n')
    f.write("w NCNO                  " + str(best_job_size) + '\n')

f.close()

# -----------------------------------------------------------------------------
# Create directories now
# -----------------------------------------------------------------------------
end_time = time.time()
#
elapsed_time = end_time - start_time
#
print("-- Single core runtime of this script:", round(elapsed_time,1), "s")
print(FMT_SEP)


