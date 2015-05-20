#!/usr/bin/python
# -*- coding: utf-8 -*-

import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs 
import numpy as np
### get parameters
execfile("./10_12s.py")

import os
if not os.path.exists(OD):
    os.makedirs(OD)
    

if not read_fields_from_file:    

    ### Readin stage


    ofr = vtk.vtkPOpenFOAMReader()
    ofr.SetFileName(ID+IF)
    ofr.SetDecomposePolyhedra(0)
    ofr.CacheMeshOn()
    ofr.SetCreateCellToPoint(0)
    ofr.DisableAllCellArrays()
    ofr.SetCellArrayStatus(fieldname,1)
    ofr.Update()

    times = dsa.vtkDataArrayToVTKArray( ofr.GetTimeValues() ,ofr)
    times = [t for t in times if t>=t0 and t<=tf]
    N = len(times)
    np.save(times_filename,times)

    cq = vtk.vtkCellQuality()
    cq.SetInputConnection(0,ofr.GetOutputPort(0))
    cq.SetQualityMeasureToVolume()
    cq.Update()
    geom = cq.GetOutputDataObject(0).GetBlock(0)

    V = np.copy(dsa.WrapDataObject(geom).CellData['CellQuality'])
    for i in range(geom.GetCellData().GetNumberOfArrays()):
        geom.GetCellData().RemoveArray(0)
    #导入POD模态数据
    Vtotal = sum(V)
    V /= Vtotal
    geom.GetCellData().AddArray(dsa.numpyTovtkDataArray(V,'vol_weight'))
    #get size of the each snapshot

    ugw = vtk.vtkXMLUnstructuredGridWriter()
    ugw.SetInputDataObject(geom)
    print 'Output to file: ',OD+'Geometry.vtu'
    ugw.SetFileName(geom_filename)
    # 默认是ascii格式，特别大！
    #ugw.SetDataModeToAscii()
    ugw.SetDataModeToBinary()
    #启用压缩
    ugw.SetCompressorTypeToZLib()
    #写
    ugw.Write()
    #disconnect cq and ofr
    cq.RemoveAllInputConnections(0)

    L = V.size
    N = len(times)
    fields = np.zeros([N,L])
    pipepout = ofr
    for i in range(N):
        t = times[i]
        print 'reading time:{}'.format(t)
        #set time value
        pipepout.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),t)
        pipepout.Update()
        d = dsa.WrapDataObject(pipepout.GetOutput().GetBlock(0))
        print 'reading field:{}'.format(fieldname)
        field = d.CellData[fieldname]
        #get the first component of composite dataset, it is the internalField
        fields[i]=np.copy(field)

    #存储数据，读取数据
    print 'write field data to file:',OD+'fields_{}.npz'.format(fieldname)
    np.savez(OD+'fields_{}.npz'.format(fieldname),fields)
else:
    fields = np.load(fields_filename)['arr_0']
    
    ugr = vtk.vtkXMLUnstructuredGridReader()
    ugr.SetFileName(geom_filename)
    ugr.Update()
    geom = ugr.GetOutputDataObject(0)
    V = np.copy(dsa.WrapDataObject(geom).CellData['vol_weight'])
    
    times = np.load(times_filename)
    assert times.shape[0] == fields.shape[0]
    assert fields.shape[1] == V.shape[0]
    N = times.shape[0] 
    L = fields.shape[1]

    print 'Read in dataset complete'
### POD section


# calculate average
field_avg = np.average(fields, axis=0)
if subtractAvg:
    fields = fields - field_avg

import modred as mr
modes, eigen_vals, eigen_vecs, correlation_mat = mr.compute_POD_matrices_snaps_method(fields.T,range(M),inner_product_weights=V,return_all=True)

if output_correlation_matrix:
    np.savetxt(POD_cm_filename,correlation_mat,delimiter=',')

if output_POD_temporal_modes: 
    print "output POD temporal modes",tm_filename
    # output temporal modes
    header_str = 'temporal modes\n'
    header_str += 'time,eigen value,singular value,normalized eigen value,accumulated normalized eigen value'
    for i in range(N):
        header_str += ',Mode{}'.format(i)
	singular_vals = eigen_vals**0.5
	POD_mode_energy_normalized = eigen_vals/correlation_mat.trace()[0,0]
    np.savetxt(POD_tm_filename, \
				np.c_[times, \
					eigen_vals**0.5, \
					eigen_vals, \
					POD_mode_energy_normalized, \
					np.cumsum(POD_mode_energy_normalized), \
					eigen_vecs], \
				delimiter = ',', \
				header = header_str)
	

if output_POD_spatial_modes:
    print "output POD spatial Modes to ",sm_filename
    #output to xml vtk unstructured grid file
    ugcd = geom.GetCellData()
    ugcd.Reset()
    ugcd.CopyAllOff()
    for i in range(ugcd.GetNumberOfArrays()):
        ugcd.RemoveArray(0)
    # import POD mode
    for i in range(M):
        ugcd.AddArray(dsa.numpyTovtkDataArray(modes[:,i],prefix+'_POD_mode_{}_{}'.format(fieldname,i)))
    # add average field
    ugcd.AddArray(dsa.numpyTovtkDataArray(field_avg,'field_{}_avg'.format(fieldname)))

    ugw = vtk.vtkXMLUnstructuredGridWriter()
    ugw.SetInputDataObject(geom)
    ugw.SetFileName(POD_sm_filename)
    ugw.Write()

print "Begin to calculate DMD modes"
modes, eigen_vals, eigen_vecs, correlation_mat = mr.compute_DMD_matrices_snaps_method(fields.T,range(M),inner_product_weights=V,return_all=True)
