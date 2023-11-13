#!/bin/python
#------------------------------------------------------------------------------
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#
# @description: 
#> Library of domain decomposition stuff
# -----------------------------------------------------------------------------
#
import numpy as np
#
# -----------------------------------------------------------------------------
# Get the section dimensions
# -----------------------------------------------------------------------------
def sec_dimensions(dims, spcng, dmn_size):
    #
    # sec_dims = np.floor(dims*spcng / dmn_size)
    #
    # -----------------------------------------------------------------------------
    # DDC uses this kind of domain decomposiiton
    # -----------------------------------------------------------------------------
    vox_per_dmn_size = np.around(dmn_size/spcng,0)
    sec_dims = np.floor(dims/vox_per_dmn_size)
    #
    return(sec_dims)
#
# -----------------------------------------------------------------------------
# Get the domain number
# -----------------------------------------------------------------------------
def get_dmn_number(dims, spcng, dmn_size, sections):
    #
    sec_dims = sec_dimensions(dims, spcng, dmn_size)
    #
    dmn_no = sections[0] + sections[1]*sec_dims[0] + sections[2]*sec_dims[0]*sec_dims[1]
    #
    return(dmn_no)
#
# -----------------------------------------------------------------------------
# Get the total number of domains in the image
# -----------------------------------------------------------------------------
def get_total_domains(dims, spcng, dmn_size):
    #
    tot_dmns = np.prod(sec_dimensions(dims, spcng, dmn_size))
    #
    return(tot_dmns)
#
# -----------------------------------------------------------------------------
# Get the section based on the domain number
# -----------------------------------------------------------------------------
def get_section(domain_no, dims, spcng, dmn_size):
    #
    section = np.zeros(3)
    #
    sec_dims = sec_dimensions(dims, spcng, dmn_size)
    #
    section[2] = np.floor(domain_no/(sec_dims[0]*sec_dims[1]))
    section[1] = np.floor((domain_no-(sec_dims[0]*sec_dims[1]*section[2]))/(sec_dims[0]))
    section[0] = np.floor(domain_no -(sec_dims[0]*sec_dims[1]*section[2]) - (sec_dims[0]*section[1]))
    #
    return(section)