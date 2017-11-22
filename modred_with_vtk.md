# modred

# windows

## 准备

安装最新的[ParaView-5.4.1-Qt5-OpenGL2-MPI-Windows-64bit.exe](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Windows&downloadFile=ParaView-5.4.1-Qt5-OpenGL2-MPI-Windows-64bit.exe)，安装位置为`C:\pv54`

安装[python2.7.3](https://www.python.org/ftp/python/2.7.3/python-2.7.3.amd64.msi)，位置`C:\Python27`

把`C:\pv54\bin\`和`C:\Python27\`加入windows的环境变量`PATH`中。

运行一下命令可以看出二者在`sys.path`上的差别：

```shell
PS P:\> python
Python 2.7.3 (default, Apr 10 2012, 23:24:47) [MSC v.1500 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> for k in sys.path: print k
...

C:\Windows\SYSTEM32\python27.zip
C:\Python27\DLLs
C:\Python27\lib
C:\Python27\lib\plat-win
C:\Python27\lib\lib-tk
C:\Python27
C:\Python27\lib\site-packages
>>> exit()
PS P:\> pvpython
Python 2.7.3 (default, Feb 13 2017, 17:21:55) [MSC v.1500 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> for k in sys.path: print k
...

C:\pv54\bin
C:\pv54\bin\python27.zip
C:\pv54\bin\DLLs
C:\pv54\bin\lib
C:\pv54\bin\lib\plat-win
C:\pv54\bin\lib\lib-tk
C:\pv54\bin
C:\pv54\bin\lib\site-packages
C:\pv54\bin\lib\site-packages\win32
C:\pv54\bin\lib\site-packages\win32\lib
>>>
```

所以需要加入环境变量`PYTHONPATH`：

```shell
C:\pv54\bin;C:\pv54\bin\DLLs;C:\pv54\bin\lib;C:\pv54\bin\lib\plat-win;C:\pv54\bin\lib\lib-tk;C:\pv54\bin;C:\pv54\bin\lib\site-packages;C:\pv54\bin\lib\site-packages\win32;C:\pv54\bin\lib\site-packages\win32\lib
```

也可以在python启动后加入：

```python
import sys
import os
#pv_root=r"C:\pv54"
pv_root=r"D:\paraview" 
pvpython_path=r"bin;bin\DLLs;bin\lib;bin\lib\plat-win;bin\lib\lib-tk;bin\lib\site-packages;bin\lib\site-packages\win32;bin\lib\site-packages\win32\lib".split(';')
for i,p in enumerate(pvpython_path):
    pvpython_path[i] = os.path.join(pv_root,p)

sys.path.extend(pvpython_path)

# import vtk works now!
```



为了安装modred，需要安装pip，可以通过下载[get-pip.py](https://bootstrap.pypa.io/get-pip.py)来完成：

```shell
PS C:\workspace\pod> python .\get-pip.py
Collecting pip
c:\users\dic17007\appdata\local\temp\tmp555h7w\pip.zip\pip\_vendor\requests\packages\urllib3\util\ssl_.py:318: SNIMissin
gWarning: An HTTPS request has been made, but the SNI (Subject Name Indication) extension to TLS is not available on thi
s platform. This may cause the server to present an incorrect TLS certificate, which can cause validation failures. You
can upgrade to a newer version of Python to solve this. For more information, see https://urllib3.readthedocs.io/en/late
st/security.html#snimissingwarning.
c:\users\dic17007\appdata\local\temp\tmp555h7w\pip.zip\pip\_vendor\requests\packages\urllib3\util\ssl_.py:122: InsecureP
latformWarning: A true SSLContext object is not available. This prevents urllib3 from configuring SSL appropriately and
may cause certain SSL connections to fail. You can upgrade to a newer version of Python to solve this. For more informat
ion, see https://urllib3.readthedocs.io/en/latest/security.html#insecureplatformwarning.
  Downloading pip-9.0.1-py2.py3-none-any.whl (1.3MB)
    100% |################################| 1.3MB 922kB/s
Collecting setuptools
  Downloading setuptools-36.7.2-py2.py3-none-any.whl (482kB)
    100% |################################| 491kB 1.9MB/s
Collecting wheel
  Downloading wheel-0.30.0-py2.py3-none-any.whl (49kB)
    100% |################################| 51kB 3.4MB/s
Installing collected packages: pip, setuptools, wheel
Successfully installed pip-9.0.1 setuptools-36.7.2 wheel-0.30.0
c:\users\dic17007\appdata\local\temp\tmp555h7w\pip.zip\pip\_vendor\requests\packages\urllib3\util\ssl_.py:122: InsecureP
latformWarning: A true SSLContext object is not available. This prevents urllib3 from configuring SSL appropriately and
may cause certain SSL connections to fail. You can upgrade to a newer version of Python to solve this. For more informat
ion, see https://urllib3.readthedocs.io/en/latest/security.html#insecureplatformwarning.
PS C:\workspace\pod>
```

然后就可以用pip安装modred了，顺带安装上numpy和future包：

```shell

PS C:\workspace\pod> python -m pip install modred
Collecting modred
C:\Python27\lib\site-packages\pip\_vendor\requests\packages\urllib3\util\ssl_.py:318: SNIMissingWarning: An HTTPS request has been made, but the SNI (Subject Name Indication) exten
sion to TLS is not available on this platform. This may cause the server to present an incorrect TLS certificate, which can cause validation failures. You can upgrade to a newer ve
rsion of Python to solve this. For more information, see https://urllib3.readthedocs.io/en/latest/security.html#snimissingwarning.
  SNIMissingWarning
C:\Python27\lib\site-packages\pip\_vendor\requests\packages\urllib3\util\ssl_.py:122: InsecurePlatformWarning: A true SSLContext object is not available. This prevents urllib3 from
 configuring SSL appropriately and may cause certain SSL connections to fail. You can upgrade to a newer version of Python to solve this. For more information, see https://urllib3.
readthedocs.io/en/latest/security.html#insecureplatformwarning.
  InsecurePlatformWarning
  Downloading modred-2.0.4.post6.tar.gz (109kB)
    100% |################################| 112kB 578kB/s
Collecting numpy (from modred)
  Downloading numpy-1.13.3-cp27-none-win_amd64.whl (13.0MB)
    100% |################################| 13.0MB 104kB/s
Collecting future (from modred)
  Downloading future-0.16.0.tar.gz (824kB)
    100% |################################| 829kB 1.2MB/s
Building wheels for collected packages: modred, future
  Running setup.py bdist_wheel for modred ... done
  Stored in directory: C:\Users\dic17007\AppData\Local\pip\Cache\wheels\bb\b5\e1\9f18da5f750a0674082c8c66651d969e60417cc941c59c83d2
  Running setup.py bdist_wheel for future ... done
  Stored in directory: C:\Users\dic17007\AppData\Local\pip\Cache\wheels\c2\50\7c\0d83b4baac4f63ff7a765bd16390d2ab43c93587fac9d6017a
Successfully built modred future
Installing collected packages: numpy, future, modred
Successfully installed future-0.16.0 modred-2.0.4.post6 numpy-1.13.3
```

测试

```python
import sys
import os
import numpy
import vtk
from vtk.numpy_interface import dataset_adapter as dsa 
import vtk.numpy_interface.algorithms as algs

data_root=r"C:\workspace\pod\cuttingPlane"
# read in the values
reader = vtk.vtkPolyDataReader()
reader.SetFileName(os.path.join(data_root,"0.01","p_zNormal.vtk"))
reader.Update()
p = dsa.WrapDataObject(reader.GetOutput()).PointData['p']

print numpy.sum(p)
## should be 2.84233e+09 
```

### vector and handle

modred推荐的使用方式是：用vector和handle包装大的计算对象。

Summarizing, to use modred on arbitrary data, define

1. A vector object that has:
   1. Vector addition (“+”, `__add__`),
   2. Scalar multiplication (“*”, `__mul__`),
   3. Optional: inherits from `vectors.Vector`.
2. A function `inner_product(vec1, vec2)`.
3. A vector handle class that has:
   1. Member function `get()` which returns a vector handle.
   2. Member function `put(vec)` where `vec` is a vector handle.
   3. Optionally inherits from `vectors.VecHandle`. If so, member function names in 1 and 2 change to `_get` and `_put`.



vector是VTKArray

handle主要是文件路径

从*.vtk文件get

put()到numpy格式的数组中。

```python
from future import standard_library
standard_library.install_hooks()
import pickle
from copy import deepcopy
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa 
import vtk.numpy_interface.algorithms as algs

## just use VTKArray as CustomVector, VTKArray is compatible with numpy.array.
import modred as mr
## just use VTKArray as CustomVector, VTKArray is compatible with numpy.array.
class CustomVector(mr.Vector):
    def __init__(self, data_array):
        self.data_array = data_array
    def __add__(self, other):
        """Return a new object that is the sum of self and other"""
        sum_vec = deepcopy(self)
        sum_vec.data_array = self.data_array + other.data_array
        return sum_vec
    def __mul__(self, scalar):
        """Return a new object that is ``self * scalar`` """
        mult_vec = deepcopy(self)
        mult_vec.data_array = mult_vec.data_array * scalar
        return mult_vec
    def inner_product(self, other):
        return np.dot(self.data_array, other.data_array)

class CustomVecHandle(mr.VecHandle):
    def __init__(self, vec_path, base_handle=None, scale=None):
        mr.VecHandle.__init__(self, base_handle, scale)
        self.vec_path = vec_path
    def _get(self):
        # read in the values
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.vec_path)
        reader.Update()
        p = dsa.WrapDataObject(reader.GetOutput()).PointData['p']
        return CustomVector(p)
    def _put(self, vec):
        np.save(self.vec_path,vec.data_array)
        return

def inner_product(v1, v2):
    return v1.inner_product(v2)
        
```



## 完整命令

下载最新的paraview 5.4.1到，安装到pv_root

安装好modred

在python中运行以下代码：

```python
pv_root=r"D:\paraview" 
data_root=r"C:\workspace\pod\cuttingPlane"
#pv_root=r"/HOME/pp569/soft/ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit" 
#data_root=r"/HOME/pp569/cheng/cuttingPlane"

## set the pythonpath to find paraview,vtk packages
import sys
import os
pvpython_path=r"bin;bin\DLLs;bin\lib;bin\lib\plat-win;bin\lib\lib-tk;bin\lib\site-packages;bin\lib\site-packages\win32;bin\lib\site-packages\win32\lib".split(';')
#pvpython_path=r"lib;lib\site-packages;lib\site-packages\vtk".split(';')
for i,p in enumerate(pvpython_path):
    pvpython_path[i] = os.path.join(pv_root,p)
    print "adding pythonpath: {}".format(pvpython_path[i])

sys.path.extend(pvpython_path)

## load vtk and paraview packages
from future import standard_library
standard_library.install_hooks()
import pickle
from copy import deepcopy
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa 
import vtk.numpy_interface.algorithms as algs
import modred as mr

## user-defined vector and handles

## common vtk object for writing
common_reader=vtk.vtkPolyDataReader()
common_reader_file = os.path.join(data_root,os.listdir(data_root)[0],"p_zNormal.vtk")
common_reader.SetFileName(common_reader_file)
common_reader.Update()
common_polydata = common_reader.GetOutput()
p_common_reader = dsa.WrapDataObject(common_polydata).PointData['p']
common_writer=vtk.vtkPolyDataWriter()
common_writer.SetInputData(common_polydata) # no SetInput in VTK 6+
common_writer_file = os.path.join(data_root,"p_test.vtk")
common_writer.SetFileName(common_writer_file)
p_common_reader[:]=1 
common_writer.Write() # return 1 means success

class CustomVector(mr.Vector):
    def __init__(self, data_array):
        self.data_array = data_array
    def __add__(self, other):
        """Return a new object that is the sum of self and other"""
        sum_vec = deepcopy(self)
        sum_vec.data_array = self.data_array + other.data_array
        return sum_vec
    def __mul__(self, scalar):
        """Return a new object that is ``self * scalar`` """
        mult_vec = deepcopy(self)
        mult_vec.data_array = mult_vec.data_array * scalar
        return mult_vec
    def inner_product(self, other):
        return np.dot(self.data_array, other.data_array)

def inner_product(v1, v2):
    return v1.inner_product(v2)


class CustomVecHandle(mr.VecHandle):
    def __init__(self, vec_path, base_handle=None, scale=None):
        mr.VecHandle.__init__(self, base_handle, scale)
        self.vec_path = vec_path
    def _get(self):
        # read in the values
        print "reading data from {}".format(self.vec_path)
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.vec_path)
        reader.Update()
        p = dsa.WrapDataObject(reader.GetOutput()).PointData['p']
        return CustomVector(p)
    def _put(self, vec):
        print "writing data to {}".format(self.vec_path)
        common_writer.SetFileName(self.vec_path)
        p_common_reader[:]=vec.data_array
        ier = common_writer.Write() # return 1 means success
        if ier!=1:
            raise Error()
        return


    
## POD data preparation
import os
dirs = [d for d in os.listdir(data_root) if  os.path.isdir(os.path.join(data_root,d))]
## eliminate non-numerical ones
dirs2 = []
for d in dirs:
    try:
        time = float(d)
    except ValueError:
        continue
    dirs2.append(d)

dirs = sorted(dirs2, key=lambda d:float(d))
times = [float(d) for d in dirs]
print "We will process those dirs:",dirs

## POD
vec_handles=[CustomVecHandle(os.path.join(data_root,d,"p_zNormal.vtk")) for d in dirs]
myPOD=mr.PODHandles(inner_product,max_vecs_per_node=20)
eigvals, eigvecs = myPOD.compute_decomp(vec_handles)
print "eigen values:"
print eigvals
print "eigen vectors:"
print eigvecs
## prepare handles for POD modes
N = 2
pod_mode_dir=os.path.join(data_root,"pod_modes")
if not os.path.exists(pod_mode_dir):
    os.makedirs(pod_mode_dir)# create dir in case of IO error

mode_handles=[CustomVecHandle(os.path.join(pod_mode_dir,"mode_{}.vtk".format(i))) for i in range(N) ]
myPOD.compute_modes(list(range(N)), mode_handles)

# projection coeffs.
proj_coeffs = myPOD.compute_proj_coeffs()
```



## 超算环境

```shell
source /HOME/pp336/scripts/set-module.sh #额外软件包
module load Python/2.7.9-gcc4.9.2
module load cmake/3.8.1
module load MPI/Gnu/MPICH/3.2-gcc4.9.2-dyn 
export CC=`which gcc` 
export CXX=`which g++`
export FC=`which gfortran` 

unzip setuptools-37.0.0.zip
cd setuptools-37.0.0
pvpython setup.py install --user
cd ..

tar xf python-future-0.16.0.tar.gz
cd python-future-0.16.0
pvpython setup.py install --user
cd ..
python -c "import future"

tar xf modred-2.0.4.post6.tar.gz
cd modred-2.0.4.post6
pvpython setup.py install --user
cd ..
python -c "import modred"


```



## 附录：Tianhe-2软件包

```shell
[pp569@ln0%tianhe2-C cuttingPlane]$ module avail

--------------------------- /BIGDATA/app/modulefiles ---------------------------
abinit/7.10.1                      mgltools/1.5.6
abinit/7.10.4                      mongodb/3.2.3
abinit/8.0.8                       mono/2.10.9
abruijn/1.0-333381f                MPB/1.5
ActiveTcl/8.6                      mpc/0.8.1
ADTEx/2.0                          MPFR/2.4.2
anaconda/1.9.2                     mpiblast/1.6.0
anaconda/2.1.0                     MPICH/Gnu/3.1
anaconda/4.2.0                     MPICH/Gnu/3.1-4.8.4
anaconda3/4.2.0                    MPICH/Gnu/3.1-4.9.2
arpack/arpack96                    MPICH/Gnu/3.2-gcc4.4.7-dyn
ARPACK/96-icc                      MPICH/Gnu/3.2-gcc4.8.5-dyn-noglex
ARPACK/96-icc15                    MPICH/Gnu/3.2-gcc4.9.2-dyn
ARWpost/3.1/00-CF-14               MPICH/Gnu/3.2-gcc5.4.0-dyn
ARWpost/3.1/01-CF-15               MPICH/Intel/2.1.5
aspect/1.3-icc15                   MPICH/Intel/3.1
atk/2.16.0                         MPICH/Intel/3.1-dbg
atlas/3.10.2                       MPICH/Intel/3.1-dyn
atompaw/4.0.0.12                   MPICH/Intel/3.1-icc11
AutoDock/4.2.6                     MPICH/Intel/3.1-icc13
AutoDock_Vina/1.1.2                MPICH/Intel/3.1-icc13-dyn
bazel/20160825                     MPICH/Intel/3.1-icc15-dyn
bcbio/0.9.9                        MPICH/Intel/3.1-large
bcftools/1.3.1                     MPICH/Intel/3.2-icc14-dyn
bcl2fastq/2.17                     MPICH/Intel/3.2-icc2017-dyn
beagle/2.1                         MPICH/Intel/3.2-icc2018-dyn
BEAST/1.8.1                        mplayer/1.3.0
BEAST/1.8.2                        MSTmap/20161226
BEAST/2.3.1                        multiz-tba/20160323
bedops/2.4.20                      muparser/2.2.5
bedtools/2.23.0                    muparser/2.2.5-icc15
bedtools/2.25.0                    mvapich2/2.2-gcc4.8.5
BerkeleyGW/1.1~beta                mvapich2/2.2-icc14
BerkeleyGW/1.2.0                   NAMD/2.10
bismark/0.16.3                     NAMD/2.10-mic
BLAS/3.5.0                         NAMD/2.11
BLAS/3.5.0-icc                     ncbi-vdb/2.6.4
BLAS/3.5.0-icc15                   ncbi-vdb/2.7.0
blast/2.2.30                       NCL/6.3.0
blast/2.5.0                        nco/4.6.0
blat/36                            nco/4.6.0-icc13
blitz/0.10                         ncview/2.1.5
boost/1.41.0                       nedit/5.6
boost/1.41.0-P279                  Nemo5/r17881
boost/1.41.0-P279-icc              Nemo5/r20038
boost/1.54.0                       Nemo5/r20175
boost/1.54.0-icc-MPI               netcdf/3.6.3
boost/1.58.0-gcc492-MPI            netcdf/4.1.3/00-CF-14
boost/1.58.0-icc-MPI               netcdf/4.1.3/01-CF-13
boost/1.59.0-gcc492-OPENMPI        netcdf/4.1.3/02-CF-14-para
boost/1.59.0-icc15-MPI             netcdf/4.1.3/03-CF-par
boost/1.62.0                       netcdf/4.3.2/01-CF-14
Bowtie/1.1.2                       netcdf/4.3.2/02-CF-13
Bowtie2/2.2.5                      netcdf/4.3.2/04-CF-4.9.2
Bowtie2/2.2.6                      netcdf/4.3.2/05-par-4.9.2
breakdancer/20160325               netcdf/4.3.3.1/01-CF-15
bwa/0.7.10                         netcdf/4.3.3.1/02-CF-15-parallel
bwa/0.7.12                         netcdf/4.4.1/00-CF-15-par
bwa/0.7.15                         ngs/1.2.4
bzip2/1.0.6                        ngs/1.2.5
cairo/1.14.2                       ngsqctoolkit/2.3.3
CAMx/6.30                          nwchem/6.5
cblas/3.5.0                        ocaml/4.02.3
cdo/1.7.1                          octave/3.8.2
cdo/1.7.1-icc13                    octopus/4.0.1
cellsys/5.0                        octopus/4.1.2
cesm/1.2.2                         octopus/5.0.1
cgnslib/3.2.1                      opam/1.2.2
Chaste/3.3                         Openbabel/2.3.2
Chaste/3.3_source                  Openbabel/2.3.2-python
chemshell/3.6.0                    OpenCASCADE/0.17
cif2cell/1.2.10                    OpenCASCADE/0.17-icc15
circexplorer2/2.2.6                opencmiss/0.3
CLAPACK/3.2.1                      opencv/1.0.0
cmake/3.0.2                        opencv/2.4.11
cmake/3.5.1                        opencv/2.4.9
CMAQ/5.0.2-cb05cl_ae5_aq           opencv/3.0.0
CNTK/1.7.2                         openexr/2.0.1
CNTK/2.0                           openexr/2.2.0
CNVnator/0.3.2                     OpenFOAM/2.1.1
code_saturne/4.0.1                 OpenFOAM/2.3.1
code_saturne/4.0.3                 OpenFOAM/2.4.0
code_saturne/4.0.3-med             OpenFOAM/2.x-SOWFA
code_saturne/4.0.4                 OpenFOAM/3.0.0
code_saturne/4.0.5                 OpenFOAM/4.0
conder/8.5.6                       openmpi/1.10.2-gcc4.9.2
CONTRA/2.0.6                       openmpi/1.10.3-gcc4.9.2
Control-FREEC/7.2                  openmpi/1.10-intel
copasi/4.15                        openmpi/1.8.3
cp2k/2.6.0                         openmpi/1.8.3-intel
cp2k/2.6.1                         openmpi/2.1.1-gcc4.9.2
cp2k/3.0                           osclean
cp2k/3.0-v2                        osenv/set2
cp2k/4.0-17116                     otf2/1.5.1
cp2k/4.1                           p4est/1.1.0
cube/4.3                           p4est/1.1.0-icc15
cufflinks/2.2.1                    p7zip/9.20.1
curl/7.49.0                        packmol/parallel
cutadapt/1.14                      packmol/serial
deal.II/8.3.0-icc15                pango/1.37.0
DL_POLY/4.07                       Panoply/4.7.0
DL_POLY/4.07-update                parafem/2.0.819
DL_POLY/4.08                       parallel/20161122
elcirc/5.4c                        ParaView/3.12.0
emacs/24.3                         ParaView/4.4.0
emboss/6.6.0                       ParaView/5.0.1-bin
esmf/6.3.0rp1-icc13-NC             ParaView/5.0.1-test
esmf/7.0.0-icc13-NC                parmetis/4.0.3-mpi
esmf/7.0.0-icc14                   parmetis/4.0.3-mpi-32bit
esmf/7.0.0-icc14-debug             pcre/8.30
esmf/7.0.0-icc14-NC                pcre/8.37
ESPResSOmd/3.2.0                   pegasus/4.6.2
ESPResSOmd/3.3.1                   perl/5.26.1
EXCAVATOR/2.2                      petsc/3.1-p8
ExomeCNVTest/0.51                  petsc/3.2-p7
fastqc/0.11.5                      petsc/3.2-p7-with-mumps
fastx_toolkit/0.0.14               petsc/3.4.4-t
FDS/6.1.2                          petsc/3.5.1-icc14
FDS/6.3.0                          petsc/3.5.4
FDTD_slution/8.11.337              petsc/3.5.4-icc14
ffmpeg/0.11.1                      petsc/3.6.2
fftw/2.1.5-double                  petsc/3.6.2-icc15
fftw/3.3.4-default                 petsc/3.6.3
fftw/3.3.4-double                  petsc/3.6.3-complex
fftw/3.3.4-double-avx              petsc/3.6.3-complex-debug
fftw/3.3.4-double-avx-sse2         petsc/3.6.3-gcc492
fftw/3.3.4-double-fma-icc15        petsc/3.6.3-icc14-4.7.4
fftw/3.3.4-gcc                     petsc/3.6.3-icc14-4.9.2
fftw/3.3.4-icc-float               petsc/3.6.3-icc15-5.3.0
fftw/3.3.4-icc-single              petsc/3.6.4-4.9.2
fftw/3.3.4-MPICH2.1.5              petsc/3.6.4-icc14-4.9.2
fftw/3.3.4-single-avx              petsc/3.7.3
fftw/3.3.4-single-avx-sse2         petsc/3.7.3-icc14
fftw/mkl-14                        petsc/3.7.5-icc14
fish/2.1.1                         phast/1.3
fish/2.2.0                         phonopy/1.11.10.4
FishingCNV/1.5.3                   phonopy/1.9.6
freesurfer/5.3.0                   phonopy/1.9.7
freetype/2.6                       PhyloCSF/20160324
fsl/5.0.9                          picard/1.129
FusionMap/20150331                 picard/2.9.0
FVCOM/3.2.1                        pixman/0.32.6
FVCOM-lib/00-icc14                 pkg-config/0.29.1/00-icc
Gamess_USA/May_01_2013_R1          plink/1.07
Gamess_USA/may2-2013               plink/1.9
GATK/3.3.0                         plumed/2.1.2
GATK/3.7                           plumed/2.2.3
gcc/4.7.4                          pnetcdf/1.6.0/00
gcc/4.8.4                          pnetcdf/1.6.1/00-icc15
gcc/4.9.2                          pnetcdf/1.6.1/01-icc13
gcc/5.2.0                          postgresql/8.4.22
gcc/5.3.0                          proj/4.9.2-icc13
gcc/5.4.0                          proj/4.9.2-icc14
gdal/2.1.0-icc14                   protobuf/20160825
gdk-pixbuf/2.31.4                  protobuf/3.1.0
geos/3.5.0-icc14                   protobuf/v160510-icc
get_homologues/20170302            PSI/4.0~b5
gflags/2.1.1                       Python/2.7.10-fPIC
gflags/2.2.1                       Python/2.7.9-fPIC
ghostscript/9.18                   Python/2.7.9-gcc4.9.2
glib/2.44.1                        Python/2.7.9-icc15
glib/2.45.3                        Python/3.4.2
glibc/2.14                         Python/3.5.1
glibc/2.14-gcc4.9                  Python/3.5.1-icc
glibc/2.15-gcc4.9                  Python/3.6.3
glibc/2.17                         QMMM/GMX
globus/6.0.1                       QMMM/LMP_QE
glog/0.3.3                         Qt/4.8.6
glog/0.3.3-gflags2.2.1             Qt/5.5.0
glue/1.46-Python2.7.9              Quantum_Espresso/5.1.1
gmap/20160404                      Quantum_Espresso/5.1.2
gmp/4.3.2                          Quantum_Espresso/5.1.2_MPI_ALL
gmt/5.2.1-icc14                    Quantum_Espresso/5.2.0_MPI
gnuplot/5.0.5                      Quantum_Espresso/5.3.0_MPI
gobject-introspection/1.44.0       Quantum_Espresso/6.0_MPI
gobject-introspection/1.46.0       R/3.1.2
grads/2.0.2                        R/3.2.1
grass/7.0.4-icc14                  R/3.2.2
gri/2.12.23/00-icc14               R/3.3.0
grib_api/1.22.0/01-icc14           R/3.3.2
gromacs/4.5.4-single               RapMap/20060321
gromacs/4.5.5-double               rar/5.3.b4
gromacs/4.5.5-single               RAxML/20161227
gromacs/5.0.4-double               redis/3.2.4
gromacs/5.0.4-mic                  rMATS/3.2.1
gromacs/5.0.4-single               root/6.03.04
gromacs/5.0.4-single-avx-256       rsem/1.2.29
gromacs/5.0.4-single-serial        rstudio/server/1.0.136
gromacs/5.0.5-double               rsync/3.0.6
gromacs/5.0.5-double-avx-256       rsync/3.1.2
gromacs/5.0.5-single               salmon/20060321
gromacs/5.1.1-single-avx-256       sambamba/0.6.3
gromacs/5.1.1-single-avx-gpu       sambamba/0.6.6
gromacs/5.1.2-single-avx-256       samtools/0.1.19
gromacs/5.1.2-single-avx-gpu       samtools/1.2
gromacs/5.1.3-single-avx-256       samtools/1.3.1
gromacs/5.1.4-single-avx-256       scala/2.11.8
gsl/1.14-icc                       scalasca/2.2
gsl/1.16                           scidb/r960
gsl/1.16-gcc                       scons/2.3.1
gsl/1.16-icc                       scons/2.5.1
gsl/2.1-icc                        scorep/1.4
gtk+/3.0.0                         sentieon/201704
hadoop/2.5.0                       serf/1.3.9
harfbuzz/0.9.41                    shapeit/v2.r790
hdf5/1.8.11/00-CF-14-parallel      siesta/4.0
hdf5/1.8.11/01-CF-15-parallel      siesta/4.0b
hdf5/1.8.12/00                     silo/4.10.2
hdf5/1.8.12/01-Gnu                 silo/4.10.2-hdf5
hdf5/1.8.12/02-icc11               silo/4.10.2-hdf5-gcc
hdf5/1.8.12/03-parallel            sionlib/1.5.4
hdf5/1.8.12/04                     slepc/3.6.3
hdf5/1.8.13/00                     SLEPc/3.6.1
hdf5/1.8.13/02-CF-14               SLEPc/3.6.1-icc15
hdf5/1.8.13/03-CF-13               snappy/1.1.3
hdf5/1.8.13/04-gcc49-parallel      snptest/2.5
hdf5/1.8.13/05-gcc49               spades/3.10.1
hdf5/1.8.13/06-CF-15               spark/1.5.1_hadoop-2.6
hdf5/1.8.13/07-gcc492-mpich3.1.3   sra-tools/2.6.4
hdf5/1.8.15/gcc492                 sra-tools/2.7.0
hdf5/1.8.17/00-CF-15-par           STAR/2.5.0
hdf5/1.8.9/00                      STAR/2.5.1
hisat2/2.0.5                       stringtie/1.3.1c
hmmer/3.1b2                        subversion/1.9.5
hypre/2.10.1-icc14                 SWIG/3.0.7
hypre/2.10.1-icc15                 SWIG/3.0.7-icc15
ilmbase/2.2.0                      szip/2.1
ImageMagick/6.9.2-5                szip/2.1-gcc49
impute2/2.3.2_static               szip/2.1-Gnu
intel-compilers/11.1               szip/2.1-icc11
intel-compilers/13.0.0             szip/2.1-icc15
intel-compilers/14.0.2             TCL/8.5.11
intel-compilers/15.0.1             TCL/8.6.3
intel-compilers/2017_update4       TCL/8.6.4
intel-compilers/2018               telemac-mascaret/v7p1
intel-compilers/mkl-14             test
intel-compilers/mkl-15             theano/0.8.2
jasper/1.900.1/00-CF-14-libpng     tinker/7.1.2
jasper/1.900.1/01-CF-15-libpng     TK/8.5.19
jasper/1.900.1/02-CF-14-shared     TK/8.6.3
java/jdk1.6.0_45                   tmux/2.0
java/jdk1.7.0_71                   tophat/2.0.14
java/jdk1.8.0_11                   trilinos/12.6.3
jpeg/6b                            trilinos/12.8.1
julia/0.4.5-icc14                  Trilinos/12.2.1-icc15
kallisto/0.42.4                    trim_galore/0.4.2
kggseq/20170318                    Trimmomatic/0.35
kggseq/20170328                    Trinity/v2.3.2
kim-api/v1.7.1                     udunits/1.12.11
lamem/20160313                     udunits/2.2.19
lamem/20170228                     udunits/2.2.19-icc13
lammps/10Aug15                     UnderWorld/1.7.0
lammps/11Sep15                     UnderWorld/1.7.0-02
lammps/15May15                     UnderWorld2/20151108
lammps/17Feb16                     UnderWorld2/20160615
lammps/17Feb16-mic                 valgrind/3.10.1
lammps/17Nov16                     VarScan/2.3.7
lammps/1Dec12                      vcftools/0.1.12b
lammps/30Sep16                     vcftools/0.1.14
lammps/30Sep16-mic                 velvet/1.2.10
lammps/31Mar17                     ViennaRNA/2.1.9
lammps/9Dec14                      ViSit/2.9.2
lammps/9Dec14-mic                  vmd/1.9.2
LAPACK/3.5.0                       vmd/1.9.2-Text-mode
LAPACK/3.5.0-icc                   VTK/5.10.1
LAPACK/3.5.0-icc15                 VTK/6.2.0
lastz/1.02.00                      VTK/7.0.0
lastz/1.03.73                      VTK/7.1.1
leveldb/1.15.0                     Wannier90/2.0.0
libevent/2.0.22                    Wannier90/2.0.1
libffi/3.2.1                       WanT/2.5.0
libgd/2.1.0                        WanT/2.5.1
libgphoto2/2.5.8                   WPS/3.6.1/00-CF-14
libint/1.1.4                       WPS/3.7.1/01-CF-15
libint/1.1.5                       WPS/3.7.1/02-CF-14
libint/2.0.3                       WPS/3.8.1/00-CF-14
libpng/1.2.56/00-icc               WPS/3.8.1/02-CF-14-esmf
libunwind/1.1                      WPS/3.9/01-CF-15
libxc/1.1.0                        WPS/3.9.0.1/01-CF-15
libxc/2.0.0                        WRF/3.6.1/00-CF-14-emreal
libxc/2.0.1                        WRF/3.6.1/01-CF-14-chem
libxc/2.2.2                        WRF/3.6.1/03-CF-14-emreal-beta-0.2
libzip/1.1.2                       WRF/3.7.1/00-CF-15-emreal
likwid/4.1.0                       WRF/3.7.1/01-CF-15-emreal-beta-0.2
lmdb/v20160511                     WRF/3.7.1/02-CF-14-emreal-beta-0.2
lzlib/1.8                          WRF/3.8.1/00-CF-14-emreal
makedepf90/2.8.8                   WRF/3.8.1/02-CF-14-esmf
mcl/14-137                         WRF/3.9/00-CF-15-emreal
mecat/1.2-3898797                  xhmm/20150320
med/3.0.8                          xxdiff/4.0.1
MEEP/1.2.1                         xz/5.2.2
Meerkat/0.185                      Yambo/3.4.1
meme/4.11.1                        yaml/0.1.5-icc14
mesa/7.9.2                         zlib/1.2.8
mesa/7.9.2-gcc                     zlib/1.2.8-icc13
METIS/5.1.0                        zlib/1.2.8-icc15
METIS/5.1.0-32bit                  zlib/1.2.8-mic
METIS/5.1.0-icc15

------------------- /HOME/pp336/modules-software/modulefiles -------------------
foam-extend/1.6                        vasp/5.4.1-impi5-O2
foam-extend/3.1(default)               vasp/5.4.1-mpi3.2-O1-scalapack
impi/4.1.3.048                         vasp/5.4.1-mpi3-O1
impi/5.0.2.044                         vasp/5.4.1-mpi3-O1-scalapack
python-tools/2.7.10(default)           vasp/5.4.1-mpi3-O1-scalapack-wannier90
vasp/5.3.5-soc-mpi3-O1                 vasp/5.4.1-mpi3-O2-scalapack-wannier90
vasp/5.3.5-soc-mpi3-O2                 vasp/5.4.1-vtst-impi4-O1-scalapack
vasp/5.4.1(default)                    vasp/5.4.1-vtst-mpi3.2-O1-scalapack
vasp/5.4.1-impi5-O1                    vasp/5.4.1-vtst-mpi3-O1-scalapack

```

