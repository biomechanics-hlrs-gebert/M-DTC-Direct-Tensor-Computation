# -----------------------------------------------------------------------------
# Makefile to build the material data process tool-chain
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      25.04.2021
# Last edit: 04.01.2022
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
trgt_vrsn="v1.0.0"
bin_name="dtc"
long_name="Direct Tensor Computation"
# -----------------------------------------------------------------------------
ifeq ($(PROVIDES_GIT),YES)
# Get git hash https://jblevins.org/log/vc
# rev = $(shell git describe --tags --always)
	rev = $(shell git rev-parse HEAD)
	trgt_vrsn = $(shell git describe --tags --abbrev=0 | tee VERSION)
else
	rev = NO_GIT_REPOSITORY
	trgt_vrsn = $(shell cat VERSION)
endif
# -----------------------------------------------------------------------------
# Check for environment
check-env:
ifeq ($(SYS_ENV),)
	@echo "-----------------------------------------------"
	@echo "-- Please source environment.sh <system> first."
	@echo "-----------------------------------------------"
else
	@echo "-----------------------------------------------"
	@echo "-- Environment to build for: "$(SYS_ENV)
	@echo "-----------------------------------------------"
	$(MAKE) all
endif
# ------------------------------------------------------------------------------

# C include paths for external libraries
inc_path_flag = -I$(METIS_INCPATH)
#
# ------------------------------------------------------------------------------

# Fortran include paths for external libraries
mod_path_flag = -I$(PETSC_INCPATH)
#
# ------------------------------------------------------------------------------

# Library paths for external libraries 
lib_path_flag = -L$(LAPACK_LIBPATH) -L$(METIS_LIBPATH) -L$(PETSC_LIBPATH)
#
# ------------------------------------------------------------------------------

# Link level extra libraries 
lll_extra = -lmetis -lpetsc -llapack -lblas -ldl
#
# ------------------------------------------------------------------------------
# Build path - choose relative or absolute paths by commenting them in/out.
# build_path = $(CURDIR)
build_path=.
export build_path
#
# ------------------------------------------------------------------------------
# Directories 
# st: "Subtree" - A git procedure to inherit another repository as some sort of
# submodule. https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202
#
# All directories are given as relative paths. Preprend "$buildpath" instead 
# of "." to change that situation.
st_path= $(build_path)/geb-lib
#
st_obj_dir = $(st_path)/obj/
st_mod_dir = $(st_path)/mod/
st_f_src_dir = $(st_path)/f-src/
#
mod_dir         = $(build_path)/mod/
obj_dir         = $(build_path)/obj/
lib_dir         = $(build_path)/lib/
linpack_src_dir = $(build_path)/linpack/
bin_dir         = $(build_path)/bin/
c_src_dir       = $(build_path)/c-src/
f_src_dir       = $(build_path)/f-src/
ext_f-src       = $(build_path)/f-src/ext-src_
#
# Directory for documentation
doc_dir  = $(build_path)/doc/
html_dir = $(build_path)/html/
tex_dir  = $(build_path)/latex/
# ------------------------------------------------------------------------------

# file extensions
mod_ext = .mod
obj_ext = .o
sho_ext = .so
c_ext   = .c
f90_ext = .f90
bin_suf = _x86_64
# ------------------------------------------------------------------------------

# Clean command 
clean_cmd = rm -rf
# ------------------------------------------------------------------------------

# Compilers
compiler   = mpif90
export compiler
c_compiler = mpicc
export c_compiler
# ------------------------------------------------------------------------------
# Programming Environment - gnu, LLVM
PE = gnu
# ------------------------------------------------------------------------------
# Compile mode - dev, prod
compile_MODE = dev
# ------------------------------------------------------------------------------
# Compile flags GNU Compiler
# The subtree structure requires two directories containing modules. 
# In this case, the program root/mod directory addressed by the -J 
# http://www.hpc.icc.ru/documentation/intel/f_ug1/fced_mod.htm
# -fbackslash does not work well with tex files (!)
ifeq ($(PE),gnu)
	f90_std_IJ     = -J$(mod_dir) -I$(st_mod_dir) $(mod_path_flag)
	f90_dev_flags  = -fdefault-integer-8 -fdefault-real-8 -finstrument-functions \
		-ggdb -o -O3 -fbacktrace -fbounds-check -Wno-conversion -Wall 
	f90_prod_flags = -fdefault-integer-8 -fdefault-real-8 -O3 -fbounds-check

	c_flags_c       = $(inc_path_flag) -finstrument-functions
	c_flags_linpack = -J$(mod_dir) -I$(st_mod_dir) -fdefault-integer-8 -g -O3 \
					  -Wall -finstrument-functions # -fopenmp

	ifeq ($(compile_MODE),prod)
		c_flags_f90 = $(f90_std_IJ) $(f90_prod_flags)
	else
		c_flags_f90 =  $(f90_std_IJ) $(f90_dev_flags)
	endif
#
# Linker flags for chain links
	link_flags = $(lib_path_flag) -L$(lib_dir) -g # -fopenmp -g # -pg 
	export glb_link_flags
endif
# ------------------------------------------------------------------------------

# Generate objects
#
c-objects =  $(obj_dir)OS$(obj_ext) \
             $(obj_dir)metis_interface$(obj_ext) 
#
# For linking and cleaning
f-objects = $(st_obj_dir)mod_global_std$(obj_ext) \
            $(st_obj_dir)mod_strings$(obj_ext) \
			$(st_obj_dir)mod_user_interaction$(obj_ext) \
			$(st_obj_dir)mod_formatted_plain$(obj_ext) \
            $(st_obj_dir)mod_math$(obj_ext) \
            $(st_obj_dir)mod_mechanical$(obj_ext) \
            $(st_obj_dir)mod_meta$(obj_ext) \
            $(obj_dir)mod_parameters$(obj_ext) \
            $(obj_dir)mod_times$(obj_ext) \
            $(obj_dir)mod_auxiliaries$(obj_ext) \
            $(obj_dir)mod_chain$(obj_ext) \
            $(obj_dir)mod_puredat$(obj_ext) \
            $(obj_dir)mod_eispack$(obj_ext) \
            $(obj_dir)mod_tensors$(obj_ext) \
            $(obj_dir)mod_metis$(obj_ext) \
            $(obj_dir)mod_OS$(obj_ext) \
            $(obj_dir)mod_linfe$(obj_ext) \
            $(obj_dir)mod_tensor_opt$(obj_ext) \
            $(obj_dir)mod_blas_1$(obj_ext) \
            $(obj_dir)mod_linpack$(obj_ext) \
            $(obj_dir)mod_mat_matrices$(obj_ext) \
            $(obj_dir)mod_decomp$(obj_ext) \
            $(obj_dir)mod_vtkio$(obj_ext) \
            $(obj_dir)mod_mesh_partitioning$(obj_ext) \
            $(obj_dir)mod_write_deck$(obj_ext) \
            $(obj_dir)mod_gen_quadmesh$(obj_ext) \
            $(obj_dir)mod_struct_preprocess$(obj_ext) \
            $(obj_dir)mod_struct_calcmat$(obj_ext) \
            $(obj_dir)mod_dtc_main_subroutines$(obj_ext) \
            $(obj_dir)struct_process$(obj_ext)
#
# For linking
pd-ld-objects = $(st_obj_dir)mod_global_std$(obj_ext) \
				$(st_obj_dir)mod_strings$(obj_ext) \
				$(st_obj_dir)mod_user_interaction$(obj_ext) \
				$(obj_dir)mod_auxiliaries$(obj_ext) \
				$(obj_dir)mod_puredat$(obj_ext) \
#
# For cleaning
pd-objects = $(obj_dir)pd_dump_leaf$(obj_ext) \
             $(obj_dir)pd_dump_tree$(obj_ext) \
             $(obj_dir)pd_leaf_diff$(obj_ext) \
             $(obj_dir)pd_leaf_to_file$(obj_ext) \
             $(obj_dir)pd_merge_branch_to_tree$(obj_ext)

# -----------------------------------------------------------------------------
# Executable
main_bin = $(bin_dir)$(bin_name)_$(trgt_vrsn)$(bin_suf)
#
# ------------------------------------------------------------------------------
# Build the st directory first
st: 
	$(MAKE) all -C $(st_path)
	@echo 
# ------------------------------------------------------------------------------
# Begin Building
all: st $(main_bin) 
#
# -----------------------------------------------------------------------------
# PureDat auxiliary executables 
pd_dump_leaf_bin            = $(bin_dir)pd_dump_leaf$(bin_suf)
pd_dump_tree_bin            = $(bin_dir)pd_dump_tree$(bin_suf)
pd_leaf_diff_bin            = $(bin_dir)pd_leaf_diff$(bin_suf)
pd_leaf_to_file_bin         = $(bin_dir)pd_leaf_to_file$(bin_suf)
pd_merge_branch_to_tree_bin = $(bin_dir)pd_merge_branch_to_tree$(bin_suf)
#
pd_aux_execs = $(pd_dump_leaf_bin) $(pd_dump_tree_bin) $(pd_leaf_diff_bin) \
               $(pd_leaf_to_file_bin) $(pd_merge_branch_to_tree_bin)
#
# -----------------------------------------------------------------------------
# Object and module dependency tree
#
# See doc/Module_and_Make_Deps.odg for a graphical overview
# -----------------------------------------------------------------------------
.PHONY: all
#
all: $(main_bin) $(pd_aux_execs) end_all
#
# -----------------------------------------------------------------------------
# C targets
$(obj_dir)%$(obj_ext):$(c_src_dir)%$(c_ext)
	@echo "----- Compiling " $< " -----"
	$(c_compiler) $(c_flags_c) -c $< -o $@
	@echo 

# -----------------------------------------------------------------------------
# Timer Module
$(obj_dir)mod_times$(obj_ext):$(f_src_dir)mod_times$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_times$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_times$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Auxiliaries Module
$(obj_dir)mod_auxiliaries$(obj_ext):$(st_mod_dir)global_std$(mod_ext) \
									$(f_src_dir)mod_auxiliaries$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_auxiliaries$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_auxiliaries$(f90_ext) -o $@
	@echo 
	
# -----------------------------------------------------------------------------
# Chain modules
$(obj_dir)mod_chain$(obj_ext):$(st_mod_dir)global_std$(mod_ext) \
							$(st_obj_dir)mod_meta$(obj_ext) \
							$(mod_dir)auxiliaries$(mod_ext) \
							$(mod_dir)timer$(mod_ext) \
							$(f_src_dir)mod_chain$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_chain$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_chain$(f90_ext) -o $@
	@echo 

$(mod_dir)chain_routines$(mod_ext):$(mod_dir)chain_variables$(mod_ext)
	$(clean_cmd) $(mod_dir)chain_routines$(mod_ext)
	@echo "----- Compiling " $(f_src_dir)mod_chain$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_chain$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# PureDat
$(obj_dir)mod_puredat$(obj_ext): $(st_mod_dir)user_interaction$(mod_ext) $(f_src_dir)mod_puredat$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_puredat$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_puredat$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# EisPack
$(obj_dir)mod_eispack$(obj_ext):$(f_src_dir)mod_eispack$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_eispack$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_eispack$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Tensor transformations
$(obj_dir)mod_tensors$(obj_ext):$(f_src_dir)mod_tensors$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_tensors$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_tensors$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Metis
$(obj_dir)mod_metis$(obj_ext):$(obj_dir)metis_interface$(obj_ext) $(f_src_dir)mod_metis$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_metis$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_metis$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# VTK-I/O
$(obj_dir)mod_vtkio$(obj_ext):$(mod_dir)auxiliaries$(mod_ext) $(f_src_dir)mod_vtkio$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_vtkio$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_vtkio$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Paramter Modules
$(obj_dir)mod_parameters$(obj_ext):$(st_mod_dir)global_std$(mod_ext) $(f_src_dir)mod_parameters$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_parameters$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_parameters$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Operating System C-call wrappers
$(obj_dir)mod_OS$(obj_ext): $(st_mod_dir)global_std$(mod_ext) \
                            $(obj_dir)OS$(obj_ext) $(f_src_dir)mod_OS$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_OS$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_OS$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Finite Element Routines
$(obj_dir)mod_linfe$(obj_ext):$(st_mod_dir)global_std$(mod_ext) \
								$(st_mod_dir)mechanical$(mod_ext) \
								$(mod_dir)puredat$(mod_ext) \
								$(f_src_dir)mod_linfe$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_linfe$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_linfe$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# LinPack Objects (Blas 1 and LinPack)
$(obj_dir)mod_blas_1$(obj_ext):$(st_mod_dir)global_std$(mod_ext) $(linpack_src_dir)mod_blas_1$(f90_ext)
	@echo "----- Compiling " $(linpack_src_dir)mod_blas_1$(f90_ext) " -----"
	$(compiler) $(c_flags_linpack) -c $(linpack_src_dir)mod_blas_1$(f90_ext) -o $@
	@echo 

$(obj_dir)mod_linpack$(obj_ext):$(st_mod_dir)global_std$(mod_ext) $(mod_dir)blas_1$(mod_ext) \
                                $(linpack_src_dir)mod_linpack$(f90_ext)
	@echo "----- Compiling " $(linpack_src_dir)mod_linpack$(f90_ext) " -----"
	$(compiler) $(c_flags_linpack) -c $(linpack_src_dir)mod_linpack$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# tensor_opt module
$(obj_dir)mod_tensor_opt$(obj_ext):$(st_mod_dir)global_std$(mod_ext)   $(mod_dir)timer$(mod_ext) \
									$(mod_dir)chain_routines$(mod_ext) $(f_src_dir)mod_tensor_opt$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_tensor_opt$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_tensor_opt$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Material matrices
$(obj_dir)mod_mat_matrices$(obj_ext):$(st_mod_dir)global_std$(mod_ext)  $(mod_dir)timer$(mod_ext) \
                                     $(mod_dir)chain_routines$(mod_ext) $(mod_dir)tensor_opt$(mod_ext) \
                                     $(mod_dir)blas_1$(mod_ext)         $(mod_dir)linpack$(mod_ext) \
                                     $(f_src_dir)mod_mat_matrices$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_mat_matrices$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_mat_matrices$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Domain decomposition module
$(obj_dir)mod_decomp$(obj_ext):$(st_mod_dir)global_std$(mod_ext)  $(mod_dir)puredat$(mod_ext) \
                               $(f_src_dir)mod_decomp$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_decomp$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_decomp$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Mesh Partitioning
$(obj_dir)mod_mesh_partitioning$(obj_ext):$(st_mod_dir)global_std$(mod_ext)   $(mod_dir)puredat$(mod_ext) \
                                          $(mod_dir)decomp$(mod_ext)          $(mod_dir)timer$(mod_ext) \
                                          $(mod_dir)chain_routines$(mod_ext)  $(mod_dir)vtkio$(mod_ext) \
                                          $(obj_dir)metis_interface$(obj_ext) $(mod_dir)metis$(mod_ext) \
                                          $(f_src_dir)mod_mesh_partitioning$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_mesh_partitioning$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_mesh_partitioning$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Deck output
$(obj_dir)mod_write_deck$(obj_ext):$(st_mod_dir)global_std$(mod_ext)   $(mod_dir)puredat$(mod_ext) \
                                   $(mod_dir)decomp$(mod_ext)          $(mod_dir)timer$(mod_ext) \
                                   $(mod_dir)chain_routines$(mod_ext)  $(mod_dir)vtkio$(mod_ext) \
                                   $(obj_dir)metis_interface$(obj_ext) $(mod_dir)metis$(mod_ext) \
                                   $(mod_dir)linfe$(mod_ext)           $(mod_dir)mesh_partitioning$(mod_ext) \
                                   $(st_mod_dir)strings$(mod_ext)      $(f_src_dir)mod_write_deck$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_write_deck$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_write_deck$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Mesh generation 
$(obj_dir)mod_gen_quadmesh$(obj_ext):$(st_mod_dir)global_std$(mod_ext)  $(mod_dir)puredat$(mod_ext) \
                                     $(mod_dir)decomp$(mod_ext)         $(mod_dir)timer$(mod_ext) \
                                     $(mod_dir)chain_routines$(mod_ext) $(f_src_dir)mod_gen_quadmesh$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_gen_quadmesh$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_gen_quadmesh$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Geometry and loadcase setup 
$(obj_dir)mod_struct_preprocess$(obj_ext):$(st_mod_dir)global_std$(mod_ext)     $(st_mod_dir)strings$(mod_ext) \
                                          $(mod_dir)decomp$(mod_ext)            $(mod_dir)timer$(mod_ext) \
                                          $(mod_dir)chain_routines$(mod_ext)    $(mod_dir)vtkio$(mod_ext) \
                                          $(obj_dir)metis_interface$(obj_ext)   $(mod_dir)metis$(mod_ext) \
                                          $(mod_dir)linfe$(mod_ext)             $(mod_dir)mesh_partitioning$(mod_ext) \
                                          $(mod_dir)puredat$(mod_ext)           $(mod_dir)gen_quadmesh$(mod_ext) \
                                          $(mod_dir)write_deck$(mod_ext)        $(f_src_dir)mod_struct_preprocess$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_struct_preprocess$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_struct_preprocess$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Calculate effective stiffness parameters 
$(obj_dir)mod_struct_calcmat$(obj_ext):$(st_mod_dir)global_std$(mod_ext)   $(st_mod_dir)formatted_plain$(mod_ext) \
										$(st_mod_dir)mechanical$(mod_ext)  $(mod_dir)tensors$(mod_ext)\
                                        $(mod_dir)puredat$(mod_ext)        $(mod_dir)timer$(mod_ext) \
                                        $(mod_dir)decomp$(mod_ext)         $(mod_dir)mat_matrices$(mod_ext) \
                                        $(mod_dir)chain_routines$(mod_ext) $(mod_dir)linfe$(mod_ext) \
                                        $(f_src_dir)mod_struct_calcmat$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_struct_calcmat$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_struct_calcmat$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Execution chain for the treatment of a single MVE
$(obj_dir)mod_dtc_main_subroutines$(obj_ext):$(st_mod_dir)global_std$(mod_ext) $(mod_dir)operating_system$(mod_ext) \
                                        $(mod_dir)puredat_com$(mod_ext)        $(mod_dir)chain_routines$(mod_ext) \
                                        $(mod_dir)linfe$(mod_ext)              $(mod_dir)gen_geometry$(mod_ext) \
                                        $(f_src_dir)mod_dtc_main_subroutines$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)mod_dtc_main_subroutines$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_dtc_main_subroutines$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# main Object 
$(obj_dir)struct_process$(obj_ext):$(st_mod_dir)global_std$(mod_ext)     $(st_mod_dir)mechanical$(mod_ext) \
                                   $(st_mod_dir)meta$(mod_ext) 			 $(st_mod_dir)meta_puredat_interface$(mod_ext) \
                                   $(st_mod_dir)strings$(mod_ext)        $(mod_dir)gen_quadmesh$(mod_ext) \
                                   $(mod_dir)auxiliaries$(mod_ext)       $(obj_dir)OS$(obj_ext) \
                                   $(mod_dir)operating_system$(mod_ext)  $(mod_dir)puredat$(mod_ext) \
                                   $(mod_dir)decomp$(mod_ext)            $(mod_dir)timer$(mod_ext) \
                                   $(mod_dir)chain_routines$(mod_ext)    $(mod_dir)vtkio$(mod_ext) \
                                   $(obj_dir)metis_interface$(obj_ext)   $(mod_dir)metis$(mod_ext) \
                                   $(mod_dir)linfe$(mod_ext)             $(mod_dir)mesh_partitioning$(mod_ext) \
                                   $(mod_dir)write_deck$(mod_ext)        $(mod_dir)gen_geometry$(mod_ext) \
                                   $(mod_dir)tensors$(mod_ext)           $(mod_dir)mat_matrices$(mod_ext) \
                                   $(mod_dir)calcmat$(mod_ext)           $(mod_dir)puredat_com$(mod_ext) \
                                   $(mod_dir)petsc_opt$(mod_ext)         $(mod_dir)dtc_main_subroutines$(mod_ext) \
								   $(f_src_dir)struct_process$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)struct_process$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)struct_process$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# PureDat auxiliary executables 
$(obj_dir)pd_dump_leaf$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext) \
                                $(f_src_dir)pd_dump_leaf$(f90_ext)
	@echo "----- Compiling " pd_dump_leaf$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_dump_leaf$(f90_ext) -o $@
	@echo 

$(obj_dir)pd_dump_tree$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext) \
                                 $(f_src_dir)pd_dump_tree$(f90_ext)
	@echo "----- Compiling " pd_dump_tree$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_dump_tree$(f90_ext) -o $@
	@echo 

$(obj_dir)pd_leaf_diff$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext) \
                                 $(f_src_dir)pd_leaf_diff$(f90_ext)
	@echo "----- Compiling " pd_leaf_diff$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_leaf_diff$(f90_ext) -o $@
	@echo 

$(obj_dir)pd_leaf_to_file$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext) \
                                 $(f_src_dir)pd_leaf_to_file$(f90_ext)
	@echo "----- Compiling " pd_leaf_to_file$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_leaf_to_file$(f90_ext) -o $@
	@echo 

$(obj_dir)pd_merge_branch_to_tree$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext) \
                                 $(f_src_dir)pd_merge_branch_to_tree$(f90_ext)
	@echo "----- Compiling " pd_merge_branch_to_tree$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_merge_branch_to_tree$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Final Link step of main 
$(main_bin): export_revision $(c-objects) $(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Final link step of $(long_name) executable'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(c-objects) $(f-objects) $(lll_extra) -o $(main_bin)
	@echo
	
# -----------------------------------------------------------------------------
# Final Link step of PureDat auxiliary executables 
$(pd_dump_leaf_bin): $(pd-ld-objects) $(obj_dir)pd_dump_leaf$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Linking PureDat auxiliary pd_dump_leaf'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(pd-ld-objects) $(obj_dir)pd_dump_leaf$(obj_ext) -o $@
	@echo 

$(pd_dump_tree_bin): $(pd-ld-objects) $(obj_dir)pd_dump_tree$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Linking PureDat auxiliary pd_dump_tree'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(pd-ld-objects) $(obj_dir)pd_dump_tree$(obj_ext) -o $@

$(pd_leaf_diff_bin): $(pd-ld-objects) $(obj_dir)pd_leaf_diff$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Linking PureDat auxiliary pd_leaf_diff'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(pd-ld-objects) $(obj_dir)pd_leaf_diff$(obj_ext) -o $@
	@echo 

$(pd_leaf_to_file_bin): $(pd-ld-objects) $(obj_dir)mod_vtkio$(obj_ext) \
	                $(obj_dir)pd_leaf_to_file$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Linking PureDat auxiliary pd_leaf_to_file'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(pd-ld-objects) \
	$(obj_dir)mod_vtkio$(obj_ext) $(obj_dir)pd_leaf_to_file$(obj_ext) -o $@
	@echo 

$(pd_merge_branch_to_tree_bin): $(pd-ld-objects) $(obj_dir)pd_merge_branch_to_tree$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Linking PureDat auxiliary pd_merge_branch_to_tree'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(link_flags) $(pd-ld-objects) $(obj_dir)pd_merge_branch_to_tree$(obj_ext) -o $@
	@echo 

# --------------------------------------------------------------------------------------------------
# Export revision
export_revision:
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Write revision and git info'
	@echo "CHARACTER(LEN=scl), PARAMETER :: longname = '$(long_name)'" > $(st_f_src_dir)include_f90/revision_meta$(f90_ext)
	@echo "CHARACTER(LEN=scl), PARAMETER :: hash = '$(rev)'" >> $(st_f_src_dir)include_f90/revision_meta$(f90_ext)
	@echo "----------------------------------------------------------------------------------"

help:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- $(long_name) make targets"
	@echo "-- Regular:  »make (all)«    - Build the $(long_name)"
	@echo "-- Cleaning: »make clean«    - Remove build files, keep the geb-lib"
	@echo "-- Cleaning: »make cleanall« - Remove all build files."
	@echo "-- Docs:     »make docs      - Build the html and the tex documentation."
	@echo "----------------------------------------------------------------------------------"

docs: 
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Beginn buiding the documentation of the $(long_name)."
	@echo "----------------------------------------------------------------------------------"
	doxygen doc/doxy.conf
	$(MAKE) pdf -C $(tex_dir)  
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Successfully build the documentation of the $(long_name)."
	@echo "----------------------------------------------------------------------------------"

cleandocs:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning html documentation"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(html_dir)/*
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning tex documentation"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(tex_dir)/*
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Documentation removed."
	@echo "----------------------------------------------------------------------------------"
	
clean:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning module directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(mod_dir)*$(mod_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning object directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(obj_dir)*$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning MAIN binary"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(main_bin)
	
cleanall: clean
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning geb-lib st"
	@echo "----------------------------------------------------------------------------------"
	$(MAKE) clean -C $(st_path)

end_all: 
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Successfully built all executables."
	@echo "----------------------------------------------------------------------------------"