# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from sys import argv


cb=float(argv[1])

solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['membrane1']
solfoam.PointArrays = ['C']
t = np.array(solfoam.TimestepValues)
N=t.size

# Properties modified on calculator1
calculator1 = pv.Calculator(Input=solfoam)
calculator1.ResultArrayName = 'cp'
calculator1.Function = 'C/%5.2e'%(cb)

# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=calculator1)

outfile=open("cpavg.dat","w")
for i in range(N):
    pv.UpdatePipeline(time=t[i], proxy=int1)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
    area    = idat.CellData['Area'].item()
    cpint   =idat.PointData['cp'].item()
    print("processing time = %e\t%e\t%e\t%e" % (t[i],area,cpint,cpint/area))
    outfile.write("%e\t%e\n"%(t[i],cpint/area))

outfile.close()
