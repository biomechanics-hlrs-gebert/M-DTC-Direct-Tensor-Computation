#!/bin/python
#------------------------------------------------------------------------------
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#
# @description: 
#> Mechanical subroutines
# -----------------------------------------------------------------------------
import math as m
import numpy as np
#
# -----------------------------------------------------------------------------
# prepare transformation
# -----------------------------------------------------------------------------
def transform(AX, AY, AZ):
    #
    # rot_m = Rotation.from_euler('xyz', [AX, AY, AZ], degrees=True).as_matrix()
    ax = AZ * m.pi / 180. # alpha
    ay = AY * m.pi / 180. # phi
    az = AX * m.pi / 180. # eta
    #
    n = np.array([ m.cos(ay)*m.sin(az), m.sin(ay)*m.sin(az), m.cos(az) ])
    #
    n = n / m.sqrt(sum(n*n))
    #
    aa = rot_alg(n, ax)
    #
    return(aa)
#
# -----------------------------------------------------------------------------
# Transform a 6x6 matrix
# -----------------------------------------------------------------------------
def transform_matrix(aa, tensor_in):
    #
    BB = tra_R6(aa)
    #
    bb = np.matrix.transpose(BB)
    #
    tensor_out = np.matmul(np.matmul(bb, tensor_in), BB)
    #
    return(tensor_out)
#
# -----------------------------------------------------------------------------
# Transform a point
# -----------------------------------------------------------------------------
def transform_point(aa, xyz_mm):
    #
    xyz_mm_rotated = np.dot(aa, xyz_mm)
    #
    return(xyz_mm_rotated)
#
# -----------------------------------------------------------------------------
# iso compliance voigt
# -----------------------------------------------------------------------------
def iso_compliance_voigt(E, v):
    t_iso_inv = np.zeros((6,6))
    fctr = 1./E

    t_iso_inv[0] = [ 1. ,    -v,    -v,         .0,         .0,         .0 ]
    t_iso_inv[1] = [  -v,   1. ,    -v,         .0,         .0,         .0 ]
    t_iso_inv[2] = [  -v,    -v,   1. ,         .0,         .0,         .0 ]
    t_iso_inv[3] = [ .0 ,    .0,    .0,  2.*(1.+v),         .0,         .0 ]
    t_iso_inv[4] = [ .0 ,    .0,    .0,         .0,  2.*(1.+v),         .0 ]
    t_iso_inv[5] = [ .0 ,    .0,    .0,         .0,         .0,  2.*(1.+v) ]

    return(t_iso_inv*fctr)
#
# -----------------------------------------------------------------------------
# Rotation matrix for R6x6
# -----------------------------------------------------------------------------
def tra_R6(aa):
    BB = np.zeros((6,6))

    BB[0] = [ aa[0, 0]**2, aa[0, 1]**2, aa[0, 2]**2, m.sqrt(2)*aa[0, 0]*aa[0, 1], m.sqrt(2)*aa[0, 0]*aa[0, 2], m.sqrt(2)*aa[0, 1]*aa[0, 2] ]
    BB[1] = [ aa[1, 0]**2, aa[1, 1]**2, aa[1, 2]**2, m.sqrt(2)*aa[1, 0]*aa[1, 1], m.sqrt(2)*aa[1, 0]*aa[1, 2], m.sqrt(2)*aa[1, 1]*aa[1, 2] ]
    BB[2] = [ aa[2, 0]**2, aa[2, 1]**2, aa[2, 2]**2, m.sqrt(2)*aa[2, 0]*aa[2, 1], m.sqrt(2)*aa[2, 0]*aa[2, 2], m.sqrt(2)*aa[2, 1]*aa[2, 2] ]

    BB[3] = [ m.sqrt(2)*aa[1, 0]*aa[0, 0], m.sqrt(2)*aa[1, 1]*aa[0, 1], m.sqrt(2)*aa[1, 2]*aa[0, 2], \
        aa[1, 0]*aa[0, 1]+aa[1, 1]*aa[0, 0], aa[1, 0]*aa[0, 2]+aa[1, 2]*aa[0, 0], aa[1, 1]*aa[0, 2]+aa[1, 2]*aa[0, 1] ]
    BB[4] = [ m.sqrt(2)*aa[0, 0]*aa[2, 0], m.sqrt(2)*aa[0, 1]*aa[2, 1], m.sqrt(2)*aa[0, 2]*aa[2, 2], \
        aa[0, 0]*aa[2, 1]+aa[0, 1]*aa[2, 0], aa[0, 0]*aa[2, 2]+aa[0, 2]*aa[2, 0], aa[0, 1]*aa[2, 2]+aa[0, 2]*aa[2, 1] ]
    BB[5] = [ m.sqrt(2)*aa[1, 0]*aa[2, 0], m.sqrt(2)*aa[1, 1]*aa[2, 1], m.sqrt(2)*aa[1, 2]*aa[2, 2], \
        aa[1, 0]*aa[2, 1]+aa[1, 1]*aa[2, 0], aa[1, 0]*aa[2, 2]+aa[1, 2]*aa[2, 0], aa[1, 1]*aa[2, 2]+aa[1, 2]*aa[2, 1] ]

    return(BB)
#
# -----------------------------------------------------------------------------
# General rotation in R3
# -----------------------------------------------------------------------------
def rot_alg(axis, angle):
    rr = np.zeros((3,3))

    rr[0,0] = m.cos(angle) + axis[0]*axis[0]* (1. - m.cos(angle))
    rr[0,1] = axis[0]*axis[1]* (1. - m.cos(angle)) - axis[2] * m.sin(angle)
    rr[0,2] = axis[0]*axis[2]* (1. - m.cos(angle)) + axis[1] * m.sin(angle)

    rr[1,0] = axis[1]*axis[0]* (1. - m.cos(angle)) + axis[2] * m.sin(angle)
    rr[1,1] = m.cos(angle) + axis[1]*axis[1]* (1. - m.cos(angle))
    rr[1,2] = axis[1]*axis[2]* (1. - m.cos(angle)) - axis[0] * m.sin(angle)

    rr[2,0] = axis[2]*axis[0]* (1. - m.cos(angle)) - axis[1] * m.sin(angle)
    rr[2,1] = axis[2]*axis[1]* (1. - m.cos(angle)) + axis[0] * m.sin(angle)
    rr[2,2] = m.cos(angle) + axis[2]*axis[2]* (1. - m.cos(angle))

    return(rr)