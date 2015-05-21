#!/usr/bin/python
# -*- coding: utf-8 -*-
#using execfile(filename) to load these settings
#####################[IOData]###########################
# working directory and data filename
ID="""/home/wei/Wu/VortexStreet/New/150cylinder13Wrerun/"""
IF="""150cylinder13W.foam"""
# output directory
OD="""/media/seagate_d4/catdog/POD_output/VortexStreet/10_12s"""
#start time
t0=12
# end time
tf=13
# name of field to be decomposed
fieldname = "p"
# prefix of spatial modes' name
prefix="10_12s_"
# True= read data from case files; False= read data from fields.npz
read_fields_from_file=True
# Variable interpolation form= %(NAME)s is available
fields_filename="""{}/fields_{}.npz""".format(OD,fieldname)
geom_filename="""{}/Geometry.vtu""".format(OD)
times_filename="""{}/times.npy""".format(OD)

do_POD = False
do_DMD = True
#####################[POD]###########################
# accuracy dete rmines how many modes will be calculated.
accuracy=0.9999
# spatial mode to be calculated, M=0 means determined by program according to accuracy
M=10
#
subtractAvg = False
#
useVolWeight = True
#
output_correlation_matrix=True
POD_cm_filename="""{}/correlation_matrix.csv""".format(OD)

output_POD_temporal_modes=True
POD_tm_filename="""{}/POD_temporal_modes.csv""".format(OD)

output_POD_spatial_modes=True
POD_sm_filename="""{}/POD_spatial_modes.vtu""".format(OD)

doReconstruction = True
#选10个点输出
MR=10 
POD_reconstruction_filename="""{}/POD_reconstruction.vtu""".format(OD)

#####################[DMD]###########################
# setting of DMD decomposition
M_DMD=10

output_DMD_info=True
DMD_info_filename="""{}/dmd_info.csv""".format(OD)

output_DMD_build_coeffs=True
DMD_build_coeffs_filename="""{}/dmd_build_coeffs.npz""".format(OD)

output_DMD_spatial_modes=True
DMD_sm_filename="""{}/DMD_spatial_modes.vtu""".format(OD)