#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Takes in fullWarp_Abs provided by grad_unwarp toolbox from Connectome
and subtracts out x/y/z coordinates to create a relative displacement
file that 3dNwarpApply can use

copied from code by cao c. March 2018
edited by mps 2020.07.01

"""
import os, nibabel, sys
import numpy as np

fullWarpAbs = sys.argv[1]
if fullWarpAbs is None:
    fullWarpAbs = input('fullWarp_abs file # (w/ single quotes): ')

f, ext = os.path.splitext(fullWarpAbs)
if ext == '.gz': # we have a .nii.gz file
    ext = '.nii.gz'

try:
    p, f = os.path.split(fullWarpAbs)
    fullWarpRel = os.path.join(p, 'fullWarp_rel' + ext)
except:
    fullWarpRel = 'fullWarp_rel' + ext
    
fullWarpPointer = nibabel.load(fullWarpAbs)
fullWarpData = fullWarpPointer.get_data()
rel = np.zeros(fullWarpData.shape)

fullWarpRAS_code = nibabel.aff2axcodes(fullWarpPointer.affine)
print('%s: (%s,%s,%s)' %(fullWarpAbs,
                         fullWarpRAS_code[0],
                         fullWarpRAS_code[1],
                         fullWarpRAS_code[2]))
                         


canonicalRAS_code = ('R', 'A', 'S') # this is nibabel canonical, they call# this is nibabel canonical, they call
# it RAS+, I would call it LPI...
# so the trick is to use fullWarpRAS_code 
# scale it to mm
voxel_size = fullWarpPointer.header.get_zooms()
if fullWarpRAS_code == ('L', 'I', 'P'): # CORONAL, RIGHT-HANDED
    x, y, z = np.meshgrid(range(fullWarpData.shape[0]-1, -1, -1), 
                          range(fullWarpData.shape[1]), 
                          range(fullWarpData.shape[2]),
                          indexing='ij')
    # Since nibabel is going to write as RAS, we need to switch y & z
    rel[:, :, :, 0] = -(fullWarpData[:, :, :, 0] - voxel_size[0]*x)
    rel[:, :, :, 1] = (fullWarpData[:, :, :, 2] - voxel_size[2]*z)
    rel[:, :, :, 2] = -(fullWarpData[:, :, :, 1] - voxel_size[1]*y)
elif fullWarpRAS_code == ('L', 'P', 'S'): # AXIAL, RIGHT-HANDED
    x, y, z = np.meshgrid(range(fullWarpData.shape[0]-1, -1, -1), 
                          range(fullWarpData.shape[1]), 
                          range(fullWarpData.shape[2]),
                          indexing='ij')
    rel[:, :, :, 0] = -(fullWarpData[:, :, :, 0] - voxel_size[0]*x)
    rel[:, :, :, 1] = fullWarpData[:, :, :, 1] - voxel_size[1]*y
    rel[:, :, :, 2] = fullWarpData[:, :, :, 2] - voxel_size[2]*z
elif fullWarpRAS_code == ('L', 'A', 'S'): # AXIAL, LEFT-HANDED
    x, y, z = np.meshgrid(range(fullWarpData.shape[0]), 
                          range(fullWarpData.shape[1]), 
                          range(fullWarpData.shape[2]),
                          indexing='ij')
    rel[:, :, :, 0] = fullWarpData[:, :, :, 0] - voxel_size[0]*x
    rel[:, :, :, 1] = -(fullWarpData[:, :, :, 1] - voxel_size[1]*y)
    rel[:, :, :, 2] = fullWarpData[:, :, :, 2] - voxel_size[2]*z
elif fullWarpRAS_code == ('L', 'S', 'P'): # CORONAL, LEFT-HANDED
    x, y, z = np.meshgrid(range(fullWarpData.shape[0]), 
                          range(fullWarpData.shape[1]), 
                          range(fullWarpData.shape[2]),
                          indexing='ij')
    # Since nibabel is going to write as RAS, we need to switch y & z
    rel[:, :, :, 0] = (fullWarpData[:, :, :, 0] - voxel_size[0]*x)
    rel[:, :, :, 1] = (fullWarpData[:, :, :, 2] - voxel_size[2]*z)
    rel[:, :, :, 2] = (fullWarpData[:, :, :, 1] - voxel_size[1]*y)
else:
    print('''!!!!!!!!!!!!!!!!!!
             We do not yet have a principled understanding of the
             interaction of slice orientation, conversion software,
             and the assumptions nibabel and 3dNwarpApply make in 
             handling these files, so each new combination needs to
             be tested empirically. The coordinate system
             detected by nibabel in fullWarp_abs header is not one
             of the cases we've worked out.
             !!!!!!!!!!!!!!!!!!
          ''')
    
# now ... when nibabel writes it, it's going to write it as RAS, regardless
# of how it came in ...
img = nibabel.Nifti1Image(rel, fullWarpPointer.affine, fullWarpPointer.header)

nibabel.save(img, fullWarpRel)
