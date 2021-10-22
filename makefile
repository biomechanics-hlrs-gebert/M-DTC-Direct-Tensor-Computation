###############################################################################
# Copyright 2021 Ralf Schneider
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################
#
#******************************************************************************
#*                                                                           **
#* makefile for the material data process tool-chain                         **
#*                                                                           **
#* last edited : on  22.10.2021                                              **
#*               by  Johannes Gebert                                         **
#*                                                                           **
#* For use of make please see                                                **
#* http://www.gnu.org/software/make/manual/make.html                         **
#*                                                                           **
#******************************************************************************
#
# -----------------------------------------------------------------------------
trgt_vrsn="v1.0.0"
bin_name="ddtc"
long_name="Directly Discretizing Tensor Computation"
# -----------------------------------------------------------------------------
# Check for environment
check-env:
ifeq ($(DDTC_ARCH),)
	@echo ----------------------------------------------------
	@echo "Please source environment.sh first."
	@echo ----------------------------------------------------
else
	@echo ----------------------------------------------------
	@echo "Environment to build for: "$(DDTC)
	$(MAKE) all
	@echo ----------------------------------------------------
	@echo "-- Successfully build $(long_name)."
	@echo ----------------------------------------------------
endif
#
#------------------------------------------------------------------------------
# C include paths for external libraries --------------------------------------
inc_path_flag = -I$(METIS_INCPATH)
#
#------------------------------------------------------------------------------
# Fortran include paths for external libraries --------------------------------
mod_path_flag = -I$(PETSC_INCPATH)
#
#------------------------------------------------------------------------------
# Library paths for external libraries ---------------------------------------- 
lib_path_flag = -L$(LAPACK_LIBPATH) -L$(METIS_LIBPATH) -L$(PETSC_LIBPATH)
#
#------------------------------------------------------------------------------
# Link level extra libraries  ------------------------------------------------- 
lll_extra = -lmetis -lpetsc -llapack -lblas -ldl
#
#------------------------------------------------------------------------------
# Build path ------------------------------------------------------------------
build_path = $(CURDIR)
export build_path
#
#------------------------------------------------------------------------------
# Directories -----------------------------------------------------------------
#
# Directory for modules -----------------------------------
mod_dir = $(build_path)/mod/
#
# Directories for objects ---------------------------------
obj_dir = $(build_path)/obj/
#
# Directories for libraries -------------------------------
lib_dir = $(build_path)/lib/
#
# Directory for binaries ----------------------------------
bin_dir = $(build_path)/bin/
#
# Directory for fortran sources ---------------------------
f_src_dir = $(build_path)/f-src/
#
# Directory for linpack sources ---------------------------
linpack_src_dir = $(build_path)/linpack/
#
# Directory for C sources ---------------------------------
c_src_dir = $(build_path)/c-src/
#
# Directory for documentation -----------------------------
doc_dir = $(build_path)/doc/
tex_dir = $(build_path)/latex/
#
#------------------------------------------------------------------------------
# Clean command ---------------------------------------------------------------
clean_cmd = rm -f
#
#------------------------------------------------------------------------------
# file extensions -------------------------------------------------------------
f90_ext  = .f90
obj_ext  = .o
mod_ext  = .mod
c_ext    =.c
#
#------------------------------------------------------------------------------
# Binary suffix ---------------------------------------------------------------
bin_suffix = _x86_64
#
#------------------------------------------------------------------------------
# Compiler, Linker and programming environment definitions --------------------
#
# Compilers -----------------------------------------------
compiler   = mpif90
export compiler
c_compiler = mpicc
export c_compiler
#
# Programming Environment ---------------------------------
PrgEnv = gnu
#
# ---------------------------------------------------------
# Compile flags GNU compiler ------------------------------
ifeq ($(PrgEnv),gnu)
   #
   # Compile flags for libraries --------------------------
   c_flags_f90     = -J$(mod_dir) -I$(mod_dir) $(mod_path_flag) \
                     -fdefault-integer-8                        \
					 -fdefault-real-8 	                        \
					 -g				                            \
				 	 -o    			                            \
	             	 -O3			                            \
	             	 -fbacktrace                                \
                 	 -fbounds-check                             \
				 	 -fbackslash                                \
	 			 	 -Wno-conversion                            \
                 	 -Wall                                      \
				 	 -finstrument-functions
                     #-fopenmp  #-finstrument-functions 
   c_flags_c       = $(inc_path_flag) \
                     -finstrument-functions
   c_flags_linpack = -J$(mod_dir) -I$(mod_dir) -fdefault-integer-8 -g -O3 \
                     -Wall \
                     -finstrument-functions # -fopenmp
   #
   # Linker flags for chain links -------------------------
   link_flags = $(lib_path_flag) -L$(lib_dir) -g # -fopenmp -g # -pg 
   export glb_link_flags
   #
endif
# ---------------------------------------------------------
# Compile flags Intel Compiler ----------------------------
ifeq ($(PrgEnv),intel)
   #
   # Compile flags for libraries --------------------------
   c_flags_f90     = -module $(mod_dir) -I$(mod_dir)  \
                     -integer-size 64 -real-size 64 -g -O3 #-fopenmp
   c_flags_c       = $(inc_path_flag)
   c_flags_linpack = -module $(mod_dir) -I$(mod_dir) -integer-size 64 -g -O3 \
                     -Wall # -fopenmp
   #
   # Linker flags for chain links -------------------------
   link_flags = $(lib_path_flag) -L$(lib_dir) -g # -fopenmp -g # -pg 
   export glb_link_flags
   #
endif
#
#------------------------------------------------------------------------------
# Objects to be generated -----------------------------------------------------
#
c-objects =  $(obj_dir)OS$(obj_ext)                    \
             $(obj_dir)metis_interface$(obj_ext) 

f-objects = $(obj_dir)mod_precision$(obj_ext)          \
            $(obj_dir)mod_parameters$(obj_ext)         \
            $(obj_dir)mod_times$(obj_ext)              \
            $(obj_dir)mod_strings$(obj_ext)            \
            $(obj_dir)mod_puredat$(obj_ext)            \
            $(obj_dir)mod_eispack$(obj_ext)            \
            $(obj_dir)mod_tensors$(obj_ext)            \
            $(obj_dir)mod_metis$(obj_ext)              \
            $(obj_dir)mod_OS$(obj_ext)                 \
            $(obj_dir)mod_linfe$(obj_ext)              \
            $(obj_dir)mod_chain$(obj_ext)              \
            $(obj_dir)mod_math$(obj_ext)               \
            $(obj_dir)mod_blas_1$(obj_ext)             \
            $(obj_dir)mod_linpack$(obj_ext)            \
            $(obj_dir)mod_mat_matrices$(obj_ext)       \
            $(obj_dir)mod_decomp$(obj_ext)             \
            $(obj_dir)mod_vtkio$(obj_ext)              \
            $(obj_dir)mod_mesh_partitioning$(obj_ext)  \
            $(obj_dir)mod_write_deck$(obj_ext)         \
            $(obj_dir)mod_gen_quadmesh$(obj_ext)       \
            $(obj_dir)mod_struct_preprocess$(obj_ext)  \
            $(obj_dir)mod_struct_calcmat$(obj_ext)     \
            $(obj_dir)struct_process$(obj_ext)

pd-objects = $(obj_dir)pd_dump_leaf$(obj_ext)            \
             $(obj_dir)pd_dump_tree$(obj_ext)            \
             $(obj_dir)pd_leaf_diff$(obj_ext)            \
             $(obj_dir)pd_leaf_to_file$(obj_ext)         \
             $(obj_dir)pd_merge_branch_to_tree$(obj_ext)

# -----------------------------------------------------------------------------
# struct-process executable ---------------------------------------------------
main_bin = $(bin_dir)struct_process$(bin_suffix)
#
# -----------------------------------------------------------------------------
# PureDat auxiliary executables -----------------------------------------------
pd_dump_leaf_bin            = $(bin_dir)pd_dump_leaf$(bin_suffix)
pd_dump_tree_bin            = $(bin_dir)pd_dump_tree$(bin_suffix)
pd_leaf_diff_bin            = $(bin_dir)pd_leaf_diff$(bin_suffix)
pd_leaf_to_file_bin         = $(bin_dir)pd_leaf_to_file$(bin_suffix)
pd_merge_branch_to_tree_bin = $(bin_dir)pd_merge_branch_to_tree$(bin_suffix)
#
pd_aux_execs = $(pd_dump_leaf_bin) $(pd_dump_tree_bin) $(pd_leaf_diff_bin) \
               $(pd_leaf_to_file_bin) $(pd_merge_branch_to_tree_bin)
#
#==============================================================================
# Object and module dependency tree
#
# See doc/Module_and_Make_Deps.odg for a graphical overview
#==============================================================================
#
.PHONY: all
#
all: $(main_bin) $(pd_aux_execs) end_all
#
# -----------------------------------------------------------------------------
#-- C targets -----------------------------------------------------------------
$(obj_dir)%$(obj_ext):$(c_src_dir)%$(c_ext)
	@echo "***** Compiling " $< " *****"
	$(c_compiler) $(c_flags_c) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Precision Module ----------------------------------------------------------
$(obj_dir)mod_precision$(obj_ext):$(f_src_dir)mod_precision$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Timer Module --------------------------------------------------------------
$(obj_dir)mod_times$(obj_ext):$(f_src_dir)mod_times$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Strings Module ------------------------------------------------------------
$(obj_dir)mod_strings$(obj_ext):$(f_src_dir)mod_strings$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- PureDat -------------------------------------------------------------------
$(obj_dir)mod_puredat$(obj_ext):$(f_src_dir)mod_puredat$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- EisPack -------------------------------------------------------------------
$(obj_dir)mod_eispack$(obj_ext):$(f_src_dir)mod_eispack$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Tensor transformations ----------------------------------------------------
$(obj_dir)mod_tensors$(obj_ext):$(f_src_dir)mod_tensors$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Metis ---------------------------------------------------------------------
$(obj_dir)mod_metis$(obj_ext):$(obj_dir)metis_interface$(obj_ext) $(f_src_dir)mod_metis$(f90_ext)
	@echo "***** Compiling " mod_metis$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_metis$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- VTK-I/O -------------------------------------------------------------------
$(obj_dir)mod_vtkio$(obj_ext):$(f_src_dir)mod_vtkio$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $< -o $@

# -----------------------------------------------------------------------------
#-- Paramter Modules ----------------------------------------------------------
$(obj_dir)mod_parameters$(obj_ext):$(mod_dir)kinds$(mod_ext) \
                                   $(f_src_dir)mod_parameters$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_parameters$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Operating System C-call wrappers ------------------------------------------
$(obj_dir)mod_OS$(obj_ext): $(mod_dir)kinds$(mod_ext) \
                            $(obj_dir)OS$(obj_ext) $(f_src_dir)mod_OS$(f90_ext)
	@echo "***** Compiling " $(f_src_dir)mod_OS$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_OS$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Finite Element Routines ---------------------------------------------------
$(obj_dir)mod_linfe$(obj_ext):$(mod_dir)kinds$(mod_ext) $(f_src_dir)mod_linfe$(f90_ext)
	@echo "***** Compiling " mod_linfe$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_linfe$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- LinPack Objects (Blas 1 and LinPack) --------------------------------------
$(obj_dir)mod_blas_1$(obj_ext):$(mod_dir)kinds$(mod_ext) $(linpack_src_dir)mod_blas_1$(f90_ext)
	@echo "***** Compiling " $< " *****"
	$(compiler) $(c_flags_linpack) -c $(linpack_src_dir)mod_blas_1$(f90_ext) -o $@

$(obj_dir)mod_linpack$(obj_ext):$(mod_dir)kinds$(mod_ext) $(mod_dir)blas_1$(mod_ext) \
                                $(linpack_src_dir)mod_linpack$(f90_ext)
	@echo "***** Compiling " mod_linpack$(f90_ext) " *****"
	$(compiler) $(c_flags_linpack) -c $(linpack_src_dir)mod_linpack$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Chain modules -------------------------------------------------------------
$(obj_dir)mod_chain$(obj_ext):$(mod_dir)kinds$(mod_ext) $(mod_dir)timer$(mod_ext) \
                              $(f_src_dir)mod_chain$(f90_ext)
	@echo "***** Compiling " mod_chain$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_chain$(f90_ext) -o $@

$(mod_dir)chain_routines$(mod_ext):$(mod_dir)chain_constants$(mod_ext) \
                                   $(mod_dir)chain_variables$(mod_ext)
	$(clean_cmd) $(mod_dir)chain_routines$(mod_ext)
	@echo "***** Compiling " mod_chain$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_chain$(f90_ext) -o $@

$(mod_dir)chain_variables$(mod_ext):$(mod_dir)chain_constants$(mod_ext)
	$(clean_cmd) $(mod_dir)chain_variables$(mod_ext)
	@echo "***** Compiling " mod_chain$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_chain$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Math module ---------------------------------------------------------------
$(obj_dir)mod_math$(obj_ext):$(mod_dir)kinds$(mod_ext) $(mod_dir)timer$(mod_ext) \
                             $(mod_dir)chain_routines$(mod_ext) $(f_src_dir)mod_math$(f90_ext)
	@echo "***** Compiling " mod_math$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_math$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Material matrices ---------------------------------------------------------
$(obj_dir)mod_mat_matrices$(obj_ext):$(mod_dir)kinds$(mod_ext)          $(mod_dir)timer$(mod_ext)   \
                                     $(mod_dir)chain_routines$(mod_ext) $(mod_dir)math$(mod_ext)    \
                                     $(mod_dir)blas_1$(mod_ext)         $(mod_dir)linpack$(mod_ext) \
                                     $(f_src_dir)mod_mat_matrices$(f90_ext)
	@echo "***** Compiling " mod_mat_matrices$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_mat_matrices$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Domain decomposition module -----------------------------------------------
$(obj_dir)mod_decomp$(obj_ext):$(mod_dir)kinds$(mod_ext) $(mod_dir)puredat$(mod_ext) \
                               $(f_src_dir)mod_decomp$(f90_ext)
	@echo "***** Compiling " mod_decomp$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_decomp$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Mesh Partitioning ---------------------------------------------------------
$(obj_dir)mod_mesh_partitioning$(obj_ext):$(mod_dir)kinds$(mod_ext)           $(mod_dir)puredat$(mod_ext) \
                                          $(mod_dir)decomp$(mod_ext)          $(mod_dir)timer$(mod_ext)   \
                                          $(mod_dir)chain_routines$(mod_ext)  $(mod_dir)vtkio$(mod_ext)   \
                                          $(obj_dir)metis_interface$(obj_ext) $(mod_dir)metis$(mod_ext)   \
                                          $(f_src_dir)mod_mesh_partitioning$(f90_ext)
	@echo "***** Compiling " mod_mesh_partitioning$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_mesh_partitioning$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Deck output ---------------------------------------------------------------
$(obj_dir)mod_write_deck$(obj_ext):$(mod_dir)kinds$(mod_ext)           $(mod_dir)puredat$(mod_ext)           \
                                   $(mod_dir)decomp$(mod_ext)          $(mod_dir)timer$(mod_ext)             \
                                   $(mod_dir)chain_routines$(mod_ext)  $(mod_dir)vtkio$(mod_ext)             \
                                   $(obj_dir)metis_interface$(obj_ext) $(mod_dir)metis$(mod_ext)             \
                                   $(mod_dir)linfe$(mod_ext)           $(mod_dir)mesh_partitioning$(mod_ext) \
                                   $(mod_dir)strings$(mod_ext)         $(f_src_dir)mod_write_deck$(f90_ext)
	@echo "***** Compiling " mod_write_deck$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_write_deck$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Mesh generation -----------------------------------------------------------
$(obj_dir)mod_gen_quadmesh$(obj_ext):$(mod_dir)kinds$(mod_ext)          $(mod_dir)puredat$(mod_ext)   \
                                     $(mod_dir)decomp$(mod_ext)         $(mod_dir)timer$(mod_ext)     \
                                     $(mod_dir)chain_routines$(mod_ext) $(f_src_dir)mod_gen_quadmesh$(f90_ext)
	@echo "***** Compiling " mod_gen_quadmesh$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_gen_quadmesh$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Geometry and loadcase setup -----------------------------------------------
$(obj_dir)mod_struct_preprocess$(obj_ext):$(mod_dir)kinds$(mod_ext)           $(mod_dir)puredat$(mod_ext)      \
                                          $(mod_dir)decomp$(mod_ext)          $(mod_dir)timer$(mod_ext)        \
                                          $(mod_dir)chain_routines$(mod_ext)  $(mod_dir)vtkio$(mod_ext)        \
                                          $(obj_dir)metis_interface$(obj_ext) $(mod_dir)metis$(mod_ext)        \
                                          $(mod_dir)linfe$(mod_ext)      $(mod_dir)mesh_partitioning$(mod_ext) \
                                          $(mod_dir)strings$(mod_ext)    $(mod_dir)gen_quadmesh$(mod_ext)      \
                                          $(mod_dir)write_deck$(mod_ext) $(f_src_dir)mod_struct_preprocess$(f90_ext)
	@echo "***** Compiling " mod_struct_preprocess$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_struct_preprocess$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- Calculate effective stiffness parameters ----------------------------------
$(obj_dir)mod_struct_calcmat$(obj_ext)::$(mod_dir)kinds$(mod_ext)          $(mod_dir)tensors$(mod_ext)       \
                                        $(mod_dir)puredat$(mod_ext)        $(mod_dir)timer$(mod_ext)         \
                                        $(mod_dir)decomp$(mod_ext)         $(mod_dir)mat_matrices$(mod_ext)  \
                                        $(mod_dir)chain_routines$(mod_ext) $(mod_dir)linfe$(mod_ext)         \
                                        $(f_src_dir)mod_struct_calcmat$(f90_ext)
	@echo "***** Compiling " mod_struct_calcmat$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)mod_struct_calcmat$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- MAIN Object ---------------------------------------------------------------
$(obj_dir)struct_process$(obj_ext):$(mod_dir)kinds$(mod_ext)            $(obj_dir)OS$(obj_ext) \
                                   $(mod_dir)operating_system$(mod_ext) $(mod_dir)puredat$(mod_ext)           \
                                   $(mod_dir)decomp$(mod_ext)           $(mod_dir)timer$(mod_ext)             \
                                   $(mod_dir)chain_routines$(mod_ext)   $(mod_dir)vtkio$(mod_ext)             \
                                   $(obj_dir)metis_interface$(obj_ext)  $(mod_dir)metis$(mod_ext)             \
                                   $(mod_dir)linfe$(mod_ext)            $(mod_dir)mesh_partitioning$(mod_ext) \
                                   $(mod_dir)strings$(mod_ext)          $(mod_dir)gen_quadmesh$(mod_ext)      \
                                   $(mod_dir)write_deck$(mod_ext)       $(mod_dir)gen_geometry$(mod_ext)      \
                                   $(mod_dir)tensors$(mod_ext)          $(mod_dir)mat_matrices$(mod_ext)      \
                                   $(mod_dir)calcmat$(mod_ext)          $(mod_dir)puredat_com$(mod_ext)       \
                                   $(mod_dir)petsc_opt$(mod_ext)         $(f_src_dir)struct_process$(f90_ext)
	@echo "***** Compiling " struct_process$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)struct_process$(f90_ext) -o $@

# -----------------------------------------------------------------------------
#-- PureDat auxiliary executables ---------------------------------------------
$(obj_dir)pd_dump_leaf$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext)   \
                                 $(f_src_dir)pd_dump_leaf$(f90_ext)
	@echo "***** Compiling " pd_dump_leaf$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_dump_leaf$(f90_ext) -o $@

$(obj_dir)pd_dump_tree$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext)   \
                                 $(f_src_dir)pd_dump_tree$(f90_ext)
	@echo "***** Compiling " pd_dump_tree$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_dump_tree$(f90_ext) -o $@

$(obj_dir)pd_leaf_diff$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext)   \
                                 $(f_src_dir)pd_leaf_diff$(f90_ext)
	@echo "***** Compiling " pd_leaf_diff$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_leaf_diff$(f90_ext) -o $@

$(obj_dir)pd_leaf_to_file$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext)   \
                                 $(f_src_dir)pd_leaf_to_file$(f90_ext)
	@echo "***** Compiling " pd_leaf_to_file$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_leaf_to_file$(f90_ext) -o $@

$(obj_dir)pd_merge_branch_to_tree$(obj_ext):$(mod_dir)puredat$(mod_ext) $(mod_dir)puredat_com$(mod_ext)   \
                                 $(f_src_dir)pd_merge_branch_to_tree$(f90_ext)
	@echo "***** Compiling " pd_merge_branch_to_tree$(f90_ext) " *****"
	$(compiler) $(c_flags_f90) -c $(f_src_dir)pd_merge_branch_to_tree$(f90_ext) -o $@

# -----------------------------------------------------------------------------
# Final Link step of MAIN -----------------------------------------------------
$(main_bin): $(c-objects) $(f-objects)
	@echo '<=============================================================================>'
	@echo '<-- Final link step of struct-process executable'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(c-objects) $(f-objects) $(lll_extra) -o $(main_bin)

# -----------------------------------------------------------------------------
# Final Link step of PureDat auxiliary executables ----------------------------
$(pd_dump_leaf_bin): $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_dump_leaf$(obj_ext)
	@echo '<=============================================================================>'
	@echo '<-- Linking PureDat auxiliary pd_dump_leaf'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_dump_leaf$(obj_ext) -o $@

$(pd_dump_tree_bin): $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_dump_tree$(obj_ext)
	@echo '<=============================================================================>'
	@echo '<-- Linking PureDat auxiliary pd_dump_tree'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_dump_tree$(obj_ext) -o $@

$(pd_leaf_diff_bin): $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_leaf_diff$(obj_ext)
	@echo '<=============================================================================>'
	@echo '<-- Linking PureDat auxiliary pd_leaf_diff'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_leaf_diff$(obj_ext) -o $@

$(pd_leaf_to_file_bin): $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)mod_vtkio$(obj_ext) \
	                $(obj_dir)pd_leaf_to_file$(obj_ext)
	@echo '<=============================================================================>'
	@echo '<-- Linking PureDat auxiliary pd_leaf_to_file'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(obj_dir)mod_puredat$(obj_ext) \
	$(obj_dir)mod_vtkio$(obj_ext) $(obj_dir)pd_leaf_to_file$(obj_ext) -o $@

$(pd_merge_branch_to_tree_bin): $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_merge_branch_to_tree$(obj_ext)
	@echo '<=============================================================================>'
	@echo '<-- Linking PureDat auxiliary pd_merge_branch_to_tree'
	@echo '<=============================================================================>'
	$(compiler) $(link_flags) $(obj_dir)mod_puredat$(obj_ext) $(obj_dir)pd_merge_branch_to_tree$(obj_ext) -o $@

# -----------------------------------------------------------------------------
# Successfull end -------------------------------------------------------------
end_all: 
	@echo '<=============================================================================>'
	@echo '<-- Successfully build all !' 
	@echo '<=============================================================================>'

# -----------------------------------------------------------------------------
# Make documentation ----------------------------------------------------------
docu:
	cd $(doc_dir); doxygen 

# -----------------------------------------------------------------------------
# Clean target chain ----------------------------------------------------------
dir_clean = \
	$(MAKE) clean -C $(link);

clean: 
	@echo '<=============================================================================>'
	@echo '<-- Cleaning module directory'
	@echo '<=============================================================================>'
	$(clean_cmd) $(mod_dir)*$(mod_ext)
	@echo '<=============================================================================>'
	@echo '<-- Cleaning C object'
	@echo '<=============================================================================>'
	$(clean_cmd) $(c-objects)
	@echo '<=============================================================================>'
	@echo '<-- Cleaning Fortran objects'
	@echo '<=============================================================================>'
	$(clean_cmd) $(f-objects)
	@echo '<=============================================================================>'
	@echo '<-- Cleaning main binary'
	@echo '<=============================================================================>'
	$(clean_cmd) $(main_bin)
	@echo '<=============================================================================>'
	@echo '<-- Cleaning aux objects'
	@echo '<=============================================================================>'
	$(clean_cmd) $(pd-objects)
	@echo '<=============================================================================>'
	@echo '<-- Cleaning aux binaries'
	@echo '<=============================================================================>'
	$(clean_cmd) $(pd_aux_execs)

clean_all: clean
	@echo '<=============================================================================>'
	@echo '<-- Cleaning documetation'
	@echo '<=============================================================================>'
	$(clean_cmd) -r $(doc_dir)/html

help:
	@echo "------------------------------------------------------------------------------"
	@echo "$(long_name) make targets"
	@echo "Regular:       »make (all)«   - Build the $(long_name)"
	@echo "Cleaning:      »make clean«   - Remove generated files, keep the config"
	@echo "Documentation: »make docs     - Build the html and the tex documentation"
	@echo "------------------------------------------------------------------------------"

docs: 
	@echo "------------------------------------------------------------------------------"
	@echo "-- Beginn buiding the documentation of the $(long_name)."
	@echo "------------------------------------------------------------------------------"
	doxygen doc/doxy.conf
	$(MAKE) pdf -C $(tex_dir) 
	@echo "------------------------------------------------------------------------------"
	@echo "-- Successfully build the documentation of the $(long_name)."
	@echo "------------------------------------------------------------------------------"
