#------------------------------------------------------------------------------
#> Analyze runtimes to BVTVs
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
import pandas as pd
import numpy as np
from scipy import stats
#
sys.path.insert(1, '/home/geb/00_bone_eval_chain/P_MOD_Python')
#
from mod_decomposition import *
from mod_mechanical import *
from mod_meta_parser import *
from mod_DTC_data_connector import *
from mod_colors import colors as c
# 
# -----------------------------------------------------------------------------
# Demo for mod_DTC_data_connector.py 
# -----------------------------------------------------------------------------
#
parser = argparse.ArgumentParser(description='Give the domain size (0.6-9.6).')
parser.add_argument('-size', action="store", dest="size", help="Input meta file).")
#
cmd_args = parser.parse_args()
#
if not bool(cmd_args.size):
    print("EE Please specify the domain size.")
    exit(1)
else: 
    try:
        size = float(cmd_args.size)
    except: 
        print("EE parsing the size failed.")
        exit(2)
#
print("-- Reading data")
#
# Header 
# Domain, Domain Size, Section x, Section y, Section z, Phys. x lo, Phys. x hi, Phys. y lo, Phys. y hi, Phys. z lo, Phys. z hi, Opt. Info, Opt. Res., Pos. alpha, Pos. eta, Pos. phi, Sym. Dev., DA, BVTV, Gray Dens., DoA Zener, DoA Gebert, mps, Spec. Norm, Micro El. Type, Macro El. Type, No. Elements, No. Nodes, t-start, t-duration, ts_InitoDmn, ts_pre-P-preallo, ts_post-P-preallo, ts_Bef-Mat-as, ts_Aft-Mat-as, ts_Bef-Solve, ts_Aft-Solve, mem_InitoDmn, mem_pre-P-preallo, mem_post-P-preallo, mem_Bef-Mat-as, mem_Aft-Mat-as, mem_Bef-Solve, mem_Aft-Solve, pR_InitoDmn, pR_pre-P-preallo, pR_post-P-preallo, pR_Bef-Mat-as, pR_Aft-Mat-as, pR_Bef-Solve, pR_Aft-Solve, size_mpi_dmn, worker_main_rank, pids_returned, S11, S21, S31, S41, S51, S61, S12, S22, S32, S42, S52, S62, S13, S23, S33, S43, S53, S63, S14, S24, S34, S44, S54, S64, S15, S25, S35, S45, S55, S65, S16, S26, S36, S46, S56, S66
if size == 0.6:
    sets_x = dtc_sets_covo("FH01-1", "06", ["000008"]) 
elif size == 1.2:
    sets_x = dtc_sets_covo("FH01-1", "12", ["000032"]) 
elif size == 2.4:
    sets_x = dtc_sets_covo("FH01-1", "24", ["000256"]) 
elif size == 4.8:
    sets_x = dtc_sets_covo("FH01-1", "48", ["000508"]) 
elif size == 7.2:
    sets_x = dtc_sets_covo("FH01-1", "72", ["004064"])
elif size == 9.6:
    sets_x = dtc_sets_covo("FH01-1", "96", ["004064"]) 

# sets_x = dtc_sets_covo("FH01-1", "96", ["004064"]) 
# sets_x = dtc_sets_covo("FH01-1", "72", ["004064"]) 
#
# Domain, PBS_JOBID, worker_main_rank, npd, ppn, ppd, runtime_preallocation, runtime_matrix_assembly, runtime_solve, runtime_total, no_nodes, no_elems, preallo, mem_InitofDomain, mem_PETScPreallo_b, mem_PETScPreallo_e, mem_assembly_b, mem_assembly_e, mem_solve_b, mem_solve_e, mem_delta_comm, mem_delta_node, mem_comm_avg, mem_comm_avg_norm, mem_node_avg, mem_node_avg_norm, pids_returned, size_mpi, no_of_compute_nodes, no_of_compute_nodes_wE, E_prepost_Ws, E_solve_Ws, E_gross_Ws, P_avg_solve_W, P_avg_gross_W, P_avg_norm_solve_W, P_avg_norm_gross_W, P_avg2_solve_W, P_avg2_gross_W, P_offset, E_ex_prepost_Ws, E_ex_solve_Ws, E_ex_gross_Ws, E_norm_prepost_Ws, E_norm_solve_Ws, E_norm_gross_Ws
if size == 0.6:
    sets_y = dtc_sets_hw("FH01-1", "06", ["000008"]) 
elif size == 1.2:
    sets_y = dtc_sets_hw("FH01-1", "12", ["000032"]) 
elif size == 2.4:
    sets_y = dtc_sets_hw("FH01-1", "24", ["000256"]) 
elif size == 4.8:
    sets_y = dtc_sets_hw("FH01-1", "48", ["000508"]) 
elif size == 7.2:
    sets_y = dtc_sets_hw("FH01-1", "72", ["004064"])
elif size == 9.6:
    sets_y = dtc_sets_hw("FH01-1", "96", ["004064"]) 
#
# -----------------------------------------------------------------------------
# Merge the files to combine the information
# -----------------------------------------------------------------------------
if size in [0.6,1.2,2.4,4.8]:
    merged = pd.merge(sets_x.data[0], sets_y.data[0], on="Domain", how='inner')

    domain   = merged["Domain"]
    bvtv     = merged["BVTV"]
    no_nodes = merged["No. Nodes"]
    no_elems = merged["No. Elements"]

    runtime_preallocation   = merged["runtime_preallocation"]
    runtime_matrix_assembly = merged["runtime_matrix_assembly"]
    runtime_solve           = merged["runtime_solve"]
    runtime_total           = merged["runtime_total"]
else:
    domain   = sets_x.data[0]["Domain"]
    bvtv     = sets_x.data[0]["BVTV"]
    no_nodes = sets_x.data[0]["No. Nodes"]

spcng = 0.0149500
total_vox = int((size/spcng)+2)**3

# for ii in range(len(no_nodes)):
#     # print(no_nodes[ii], runtime_total[ii], no_nodes[ii]/runtime_total[ii])
#     try:
#         print(int(no_nodes[ii]/bvtv[ii]))
#     except:
#         continue

print("-- Voxels of monolithic domains:", total_vox)

nds_per_vox_avg = 0
nds_per_vox = 0
for ii in range(len(no_nodes)):
    # print(no_nodes[ii], runtime_total[ii], no_nodes[ii]/runtime_total[ii])
    try:
        nds_per_vox = no_nodes[ii] / (total_vox*bvtv[ii])

        if math. isinf(nds_per_vox):
            continue

        nds_per_vox_avg += nds_per_vox
        # print("-- Domain", domain[ii], " - ", round(nds_per_vox,3), "Nodes per Voxel")
    except:
        continue
print("\n\n\n")
print("-- Average nodes per Voxel: ", nds_per_vox_avg / len(no_nodes))


slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,bvtv)
print("-- r_value no_nodes,bvtv: ", r_value)

if size in [0.6,1.2,2.4,4.8]:

    slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,no_elems)
    print("-- r_value no_nodes,no_elems: ", r_value)

    slope, intercept, r_value, p_value, std_err = stats.linregress(bvtv,no_elems)
    print("-- r_value bvtv,no_elems: ", r_value)

    slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,runtime_preallocation)
    print("-- r_value no_nodes,runtime_preallocation: ", r_value)

    slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,runtime_matrix_assembly)
    print("-- r_value no_nodes,runtime_matrix_assembly: ", r_value)

    slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,runtime_solve)
    print("-- r_value no_nodes,runtime_solve: ", r_value)

    slope, intercept, r_value, p_value, std_err = stats.linregress(no_nodes,runtime_total)
    print("-- r_value no_nodes,runtime_total: ", r_value)
