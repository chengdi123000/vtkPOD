#Todo#

1.  separate code to stages, each stage should save some files for reuse
    *   Readin stage, output *fields.npz* and *grid.vtu*
    *   POD stage, output *correlation_matrix.csv*, *POD_temporal_modes.csv* and *POD_spatial_modes.vtu*
    *   DMD stage, output *DMD_eigen_values.csv* and *DMD_spatial_modes.vtu*
2.  write a configuration file to control each stage, using ConfigParser module in python to read and parse it
3.  add slice and clip algorithm to reduce size of dataset

