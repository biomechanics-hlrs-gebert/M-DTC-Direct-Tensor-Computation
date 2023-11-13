# -----------------------------------------------------------------------------
# Makefile to build the material data process tool-chain
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      17.09.2023
# Last edit: 17.09.2023
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
bin_name="create_test_case"
long_name="Direct Tensor Computation | Testcases"
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
st_path= $(build_path)/../geb-lib
#
st_obj_dir = $(st_path)/obj/
st_mod_dir = $(st_path)/mod/
st_f_src_dir = $(st_path)/f-src/
#
mod_dir   = $(build_path)/mod/
obj_dir   = $(build_path)/obj/
bin_dir   = $(build_path)/bin/
f_src_dir = $(build_path)/f-src/
#
# ------------------------------------------------------------------------------
# file extensions
# -----------------------------------------------------------------------------
mpi_ext = .mpi
mod_ext = .mod
obj_ext = .o
sho_ext = .so
f90_ext = .f90
bin_suf = _x86_64
#
# ------------------------------------------------------------------------------
# Clean command 
clean_cmd = rm -rf
# ------------------------------------------------------------------------------
#
# Inlcude build options
include $(st_path)/make.flags
#
# -----------------------------------------------------------------------------
# Fortran objects
# -----------------------------------------------------------------------------
f-objects = $(st_obj_dir)mod_global_std$(obj_ext) \
			$(st_obj_dir)mod_strings$(obj_ext) \
			$(st_obj_dir)mod_math$(obj_ext) \
			$(st_obj_dir)mod_user_interaction$(obj_ext) \
			$(st_obj_dir)mod_ser_binary$(obj_ext)\
			$(st_obj_dir)mod_formatted_plain$(obj_ext) \
			$(obj_dir)create_testcase$(obj_ext)
#
# -----------------------------------------------------------------------------
# Bionary objects
# -----------------------------------------------------------------------------
main_bin = $(bin_dir)$(bin_name)$(bin_suf)
#
# ------------------------------------------------------------------------------
# Build the st directory first
# -----------------------------------------------------------------------------
st: 
	$(MAKE) all -C $(st_path)
	@echo 
#
# ------------------------------------------------------------------------------
# Begin Building
# -----------------------------------------------------------------------------
all: st $(main_bin) 
#
# -----------------------------------------------------------------------------
# main Object 
# -----------------------------------------------------------------------------
$(obj_dir)create_testcase$(obj_ext):$(st_mod_dir)global_std$(mod_ext) \
								   $(st_mod_dir)ser_binary$(mod_ext)\
								   $(f_src_dir)create_testcase$(f90_ext)
	@echo "----- Compiling " $(f_src_dir)create_testcase$(f90_ext) "-----"
	$(f90_compiler) $(c_flags_f90)  -c $(f_src_dir)create_testcase$(f90_ext) -o $@
	@echo 
#
# -----------------------------------------------------------------------------
# Final Link step of main 
# -----------------------------------------------------------------------------
$(main_bin): $(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Final link step of $(long_name) executable'
	@echo "----------------------------------------------------------------------------------"
	$(f90_compiler) $(f-objects) -o $(main_bin)
	@echo

#
# -----------------------------------------------------------------------------
# Remove stuff
# -----------------------------------------------------------------------------
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
