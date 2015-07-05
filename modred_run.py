#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
#load parameters
execfile(sys.argv[1])

def main(argv):
    import vtk
    from vtk.numpy_interface import dataset_adapter as dsa
    from vtk.numpy_interface import algorithms as algs 
    import numpy as np
    ### get parameters
    
    import os
    if not os.path.exists(OD):
        os.makedirs(OD)
    print '!!!Output to DIR: ',OD    

    if not read_fields_from_file:    

        ### Readin stage
        # using parallel openfoam reader
        ofr = vtk.vtkPOpenFOAMReader()
        # set reader's options 
        ofr.SetFileName(ID+IF)
	print '!!!open file: ',ID+IF
        ofr.SetDecomposePolyhedra(0)
        ofr.CacheMeshOn()
        ofr.SetCreateCellToPoint(0)
        ofr.DisableAllCellArrays()
        ofr.SetCellArrayStatus(fieldname,1)
        ofr.Update()

        # VTKArray is same as numpy array
        times = dsa.vtkDataArrayToVTKArray( ofr.GetTimeValues() ,ofr)
        # select the timestep between t0 and tf
        times = [t for t in times if t>=t0 and t<=tf]
	print '!!!available time steps: ',times
        N = len(times)
        np.save(times_filename,times)

        # using CellQuality to get cell's volumes as weight
        cq = vtk.vtkCellQuality()
        cq.SetInputConnection(0,ofr.GetOutputPort(0))
        cq.SetQualityMeasureToVolume()
        cq.Update()
        # cq is a composite dataset so I need GetBlock(0)
        geom = cq.GetOutputDataObject(0).GetBlock(0)

        # get volumes of cells V, size = L (number of cells)
        V = np.copy(dsa.WrapDataObject(geom).CellData['CellQuality'])
        # normalize it as weight 
        Vtotal = sum(V)
        V /= Vtotal
        
        # delete all other CellDataArray in geom DataSet, preserve its mesh structure and topology structure
        for i in range(geom.GetCellData().GetNumberOfArrays()):
            geom.GetCellData().RemoveArray(0)
        # add volume weight to it for saving
        geom.GetCellData().AddArray(dsa.numpyTovtkDataArray(V,'vol_weight'))

        # using *.vtu file format to save the vol_weight
        ugw = vtk.vtkXMLUnstructuredGridWriter()
        ugw.SetInputDataObject(geom)
        print '!!!Output vol_weight to file: ',geom_filename
        ugw.SetFileName(geom_filename)
        # using binary format
        ugw.SetDataModeToBinary()
        # enable compression
        ugw.SetCompressorTypeToZLib()
        # write to the file
        ugw.Write()
        # disconnect cq and ofr in order to isolate this dataset object from Update()
        cq.RemoveAllInputConnections(0)

        L = V.size # number of cells
        N = len(times) #number of timesteps
        # vector data is larger in size
        if field_is_vector == True:
            fields = np.zeros([N,L,3])
        else:
            fields = np.zeros([N,L])
        pipepout = ofr
        for i in range(N):
            t = times[i]
            print '!!!reading time:{}'.format(t)
            # set time value
            pipepout.GetOutputInformation(0).Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(),t)
            # read in field data of new timestep
            pipepout.Update()
            # 
            d = dsa.WrapDataObject(pipepout.GetOutput().GetBlock(0))
            print '!!!reading field:{}'.format(fieldname)
            field = d.CellData[fieldname]
            # get the first component of composite dataset, it is the internalField
            fields[i]=np.copy(field)

        # write data to file
        print '!!!write field data to file:',fields_filename
        np.savez(fields_filename,fields)
    else: #read fields from file
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
    
    if do_POD:
        # if field is a vector, reshape the fields and corresponding volument weight
        if field_is_vector:
            shp_vec = fields.shape
            shp_flat = (fields.shape[0],fields.shape[1]*fields.shape[2])
            fields = fields.reshape(shp_flat)
            V = np.tile(V,shp_vec[2])

        # POD
        print '!!!Doing POD analysis'
        modes, eigen_vals, eigen_vecs, correlation_mat = mr.compute_POD_matrices_snaps_method(fields.T,range(M),inner_product_weights=V,return_all=True)

        # if field is a vector, reshape the output matrix
        if field_is_vector:
            fields = fields.reshape(shp_vec)
            modes = np.asarray(modes).T.reshape((modes.shape[1],shp_vec[1],shp_vec[2]))
            V = V[:shp_vec[1]]

        if output_correlation_matrix:
            print "!!!output POD correlation matrix",POD_cm_filename
            np.savetxt(POD_cm_filename,correlation_mat,delimiter=',')

        if output_POD_temporal_modes: 
            print "!!!output POD temporal modes",POD_tm_filename
            # output temporal modes
            singular_vals = eigen_vals**0.5
            POD_mode_energy_normalized = eigen_vals/correlation_mat.trace()[0,0]
            cumsum_POD_mode_energy_normalized = np.cumsum(POD_mode_energy_normalized)
            # generate header string
            header_str = 'temporal modes\n'
            header_str += 'time,eigen value,singular value,normalized eigen value,accumulated normalized eigen value'
            for i in range(N-1):
                header_str += ',Mode{}'.format(i)
            header_str += '\n'
            for i in range(N-1):
                header_str += ',SV ={}'.format(singular_vals[i])
            header_str += '\n'
            for i in range(N-1):
                header_str += ',EV ={}'.format(eigen_vals[i])
            header_str += '\n'
            for i in range(N-1):
                header_str += ',NEnergy ={}'.format(POD_mode_energy_normalized[i])
            header_str += '\n'
            for i in range(N-1):
                header_str += ',CumsumEnergy ={}'.format(cumsum_POD_mode_energy_normalized[i])
            header_str += '\n'

            np.savetxt(POD_tm_filename, \
                        np.c_[times, \
                            eigen_vecs], \
                        delimiter = ',', \
                        header = header_str)
            

        if output_POD_spatial_modes:
            print "!!!output POD spatial Modes to ",POD_sm_filename
            #output to xml vtk unstructured grid file
            ugcd = geom.GetCellData()
            ugcd.Reset()
            ugcd.CopyAllOff()
            for i in range(ugcd.GetNumberOfArrays()):
                ugcd.RemoveArray(0)
            # import POD mode
            for i in range(M):
                ugcd.AddArray(dsa.numpyTovtkDataArray(modes[i],prefix+'_POD_mode_{}_{}'.format(fieldname,i)))
            # add average field
            ugcd.AddArray(dsa.numpyTovtkDataArray(field_avg,'field_{}_avg'.format(fieldname)))

            ugw = vtk.vtkXMLUnstructuredGridWriter()
            ugw.SetInputDataObject(geom)
            ugw.SetFileName(POD_sm_filename)
            ugw.Write()
        if doReconstruction:
            print "!!! do Reconstrution with {} POD modes at time {}".format(MR,ReconTime)
            #get an empty mesh
            ugcd = geom.GetCellData()
            ugcd.Reset()
            ugcd.CopyAllOff()
            for i in range(ugcd.GetNumberOfArrays()):
                ugcd.RemoveArray(0)
            # reconstruct from first MR POD modes
            # 
            ReconN = np.searchsorted(times,ReconTime)
            print "!!!actually, reconstruction is done at time {} rather than time {}".format(times[ReconN],ReconTime)
            recon_field = np.einsum("i...,i,i",modes[:MR],eigen_vals[:MR]**0.5,np.asarray(eigen_vecs)[ReconN,:MR])+field_avg;
            ugcd.AddArray(dsa.numpyTovtkDataArray(recon_field,prefix+'_POD_{}_Reconstructed_{}_{}'.format(MR,fieldname,ReconTime)))

            ugw = vtk.vtkXMLUnstructuredGridWriter()
            ugw.SetInputDataObject(geom)
            ugw.SetFileName(POD_reconstruction_filename)
            ugw.Write()
    if do_DMD:
        print "!!!Begin to calculate DMD modes"
        # if field is a vector, reshape the fields and corresponding volument weight
        if field_is_vector:
            shp_vec = fields.shape
            shp_flat = (fields.shape[0],fields.shape[1]*fields.shape[2])
            fields = fields.reshape(shp_flat)
            V = np.tile(V,shp_vec[2])

        # DMD, I do not know which mode is important, so I have to discard modes_
        modes_, ritz_vals, mode_norms, build_coeffs = mr.compute_DMD_matrices_snaps_method(fields.T,[],inner_product_weights=V,return_all=True)

        # if field is a vector, reshape the fields, V and output matrix
        if field_is_vector:
            fields = fields.reshape(shp_vec)
            V = V[:shp_vec[1]]
        # sorting
        eorder = np.argsort(mode_norms)[::-1]
        # re-order the outputs
        ritz_vals = ritz_vals[eorder]
        mode_norms = mode_norms[eorder]
        build_coeffs = build_coeffs[:,eorder]
        #build the DMD_modes
        DMD_modes = np.einsum('ijk,il->ljk', fields,build_coeffs[:,:M_DMD])
        
        if output_DMD_info:
            print "!!!output DMD info to :",DMD_info_filename
            # output modes info
            header_str = 'DMD modes info\n'
            header_str += 'ritz_vals.real,ritz_vals.imag,growth_rate, frequency, mode_norms\n'
            header_str += r'AU,AU,1/s, Hz, AU'
            dt = np.average(times[1:]-times[:-1]) #time step
            np.savetxt(DMD_info_filename, \
                        np.c_[ np.real(ritz_vals), \
                            np.imag(ritz_vals), \
                            np.log(np.abs(ritz_vals))/dt, \
                            np.angle(ritz_vals)/dt, \
                            mode_norms], \
                        delimiter = ',', \
                        header = header_str) 

        if output_DMD_build_coeffs:
            print "!!!output DMD build coeffs. to :",DMD_build_coeffs_filename
            np.savez(DMD_build_coeffs_filename, build_coeffs)
            
        if output_DMD_spatial_modes:
            print "!!!output DMD info to :",DMD_sm_filename
            #output to xml vtk unstructured grid file
            ugcd = geom.GetCellData()
            ugcd.Reset()
            ugcd.CopyAllOff()
            for i in range(ugcd.GetNumberOfArrays()):
                ugcd.RemoveArray(0)
            #import pi
            from numpy import pi
            
            for i in range(M_DMD):
                ugcd.AddArray(dsa.numpyTovtkDataArray(np.abs(DMD_modes[i]),prefix+'_DMD_mode_abs_{}_{}'.format(fieldname,i)))
                ugcd.AddArray(dsa.numpyTovtkDataArray(np.angle(DMD_modes[i])*180/pi,prefix+'_DMD_mode_angle_{}_{}'.format(fieldname,i)))


            ugw = vtk.vtkXMLUnstructuredGridWriter()
            ugw.SetInputDataObject(geom)
            ugw.SetFileName(DMD_sm_filename)
            ugw.Write()
        
if __name__ == "__main__":
    import sys
    main(sys.argv)
