#!/usr/bin/python
# -*- coding: utf-8 -*-

import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs 
import numpy as np
import os

#start time
t0 = 10
# end time
tf = 11.026
# working directory
wd = '/home/wei/Wu/VortexStreet/New/150cylinder13Wrerun'
# filename
filename = '150cylinder13W.foam'
# name of field to be decomposed
fieldname = 'p'
# accuracy determines how many modes will be calculated.
accuracy = 0.9999
# prefix of all output file of generated data
prefix = 'volume_weighted_no_avg_substraction_'
# some switch
# True: read data from case files; False: read data from fields.npz
readField = True
#
subtractAvg = False
#
useVolWeight = True


os.chdir(wd)
print 'change to working directory:',wd
ofr = vtk.vtkPOpenFOAMReader()
ofr.SetFileName(filename)
print 'reading case: ',filename
#Set and get case type. 0 = decomposed case, 1 = reconstructed case.
ofr.SetCaseType(1) 
#for reading point/face/cell-Zones
ofr.ReadZonesOff()
#不将多面体分割成四面体
ofr.SetDecomposePolyhedra(0)
# mesh会被cache
ofr.CacheMeshOn()
#不产生点数据
ofr.SetCreateCellToPoint(0)
#只读取fieldname对应的数据
ofr.DisableAllCellArrays()
ofr.SetCellArrayStatus(fieldname,1)
print 'reading field: ',fieldname
#更新数据
ofr.Update()

#取得时间数据
times = dsa.vtkDataArrayToVTKArray( ofr.GetTimeValues() ,ofr)
#选择时间数据
times = [t for t in times if t>=t0 and t<=tf]
N = len(times)
print 'number of timesteps: ',N
print 'timesteps:',times
# test if times are evenly spaced.
if std(times[1:]-times[:-1])/average(times[1:]-times[:-1])> 1e-5:
    raise ValueError( 'The times are not evenly spaced!\nquit!')
    exit()
    
#利用vtkCellQuality Filter来计算单元体积
cq = vtk.vtkCellQuality()
#cq承接从ofr传来的输出
cq.SetInputConnection(0,ofr.GetOutputPort(0))
#设置为计算体积
cq.SetQualityMeasureToVolume()
#计算
cq.Update()
#取得计算的体积，硬拷贝一份
V = np.copy(dsa.WrapDataObject(cq.GetOutputDataObject(0)).CellData['CellQuality'].Arrays[0])
#求和取得总体积
Vtotal = sum(V)
V /= Vtotal #归一化，作为权重
L = V.size

#断开cq与ofr的Pipeline
cq.RemoveAllInputConnections(0)

if readField:
    print 'read all fields from files'
    #取得所有场数据
    fields = np.zeros([N,L])
    for i in range(N):
        t = times[i]
        print 'reading time:{}'.format(t) 
        #设置读取数据的时间
        #ofr.SetTimeValue(t)
        ofr.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),t)
        ofr.Update()
        d = dsa.WrapDataObject(ofr.GetOutput())
        print 'reading field:{}'.format(fieldname)
        field = d.CellData[fieldname]
        #get the first component of composite dataset, it is the internalField
        fields[i]=np.copy(field.Arrays[0])
        assert np.any(np.isnan(fields[i])) == False

    #存一个数据出来看看
    #ofr.GetOutputDataObject(0).GetBlock(0).GetCellData().GetArray(0).GetValue(0)    
    #设置读取数据的时间
    #ofr.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),times[-5])

    #存储数据，启用压缩
    np.savez('fields_{}.npz'.format(fieldname),fields)
else:
    #载入数据
    print 'load field from fields.npz'
    fields = np.load('fields.npz')

#减去平均值
field_avg = np.average(fields, axis=0)
if subtractAvg:
    fields = fields - field_avg


#采用稀疏线性代数方法直接求解最大的几个模态，不过似乎速度不快，算出来也不对
# 不进行空间加权
# from scipy.sparse.linalg import svds
# Etot = np.norm(fields) 
# M = max(10,N/10)
# u,s,v = svds(fields,M) 

from itertools import combinations_with_replacement as cwr

C = np.zeros([N,N])
for i,j in cwr(range(N),2):
    if useVolWeight:
        C[i,j] = np.einsum('i,i,i',fields[i],fields[j],V)
    else:
        C[i,j] = np.dot(fields[i],fields[j])
    print "correlation coeff.({},{}) calculated".format(i,j)
    if i != j:
        C[j,i] = C[i,j]

# 采用对称矩阵的特征值求解器
from numpy.linalg import eigh
ev,eVec = eigh(C)
# 有时候会有小于0的小量
ev[ev<0.0]=0.0
#ev[i]对应eVec[i]
# transpose
eVec = eVec.T

# 排序
eorder = np.argsort(ev)[::-1]
# 从大到小
ev = ev[eorder]
eVec = eVec[eorder]

sv = np.sqrt(ev)

# 计算累计分布
cen = np.cumsum(ev)/C.trace()

# 输出
for i in range(N):
    print 'Mode{},ev={},cen={}'.format(i,ev[i],cen[i])
    
# 输出时间模态和特征值
header_str = 'temporal modes\n'
header_str += 'time,eigen value,singular value,accumulated normalized eigen value'
for i in range(N):
    header_str += ',Mode{}'.format(i)
np.savetxt(prefix+'temporal_modes.csv',np.c_[times,ev,sv,cen,eVec.T],delimiter = ',',header = header_str)

# 计算达到accuracy需要的模态数
from bisect import bisect_left
M = bisect_left(cen,accuracy)
print 'Spatial Mode number to be calculated: ',M

# 计算空间模态
Modes = np.zeros([M,L])
for i in range(M):
    print 'calculating spatial mode {}'.format(i)
    Modes[i] = np.einsum( 'j,jk',eVec[i]/sv[i],fields ) 

#output to xml vtk unstructured grid file
ug = ofr.GetOutputDataObject(0).GetBlock(0)
ugcd = ug.GetCellData()
ugcd.Reset()
ugcd.CopyAllOff()
#把原来的数据删除
for i in range(ugcd.GetNumberOfArrays()):
    ugcd.RemoveArray(0)
#导入POD模态数据
for i in range(M):
    ugcd.AddArray(dsa.numpyTovtkDataArray(Modes[i],'POD_mode_{}_{}'.format(fieldname,i)))
#把平均场也加进去
ugcd.AddArray(dsa.numpyTovtkDataArray(field_avg,'field_{}_avg'.format(fieldname)))
    
# 用writer写出来
ugw = vtk.vtkXMLUnstructuredGridWriter()
ugw.SetInputDataObject(ug)
#文件名
spatial_mode_output_filename = prefix+'POD_{}_spatial_modes.vtu'.format(fieldname)
ugw.SetFileName(spatial_mode_output_filename)
# 默认是ascii格式，特别大！
#ugw.SetDataModeToAscii()
ugw.SetDataModeToBinary()
#启用压缩
ugw.SetCompressorTypeToZLib()
#写
ugw.Write()
print 'output to file:{}'.format(spatial_mode_output_filename)
print 'END'

