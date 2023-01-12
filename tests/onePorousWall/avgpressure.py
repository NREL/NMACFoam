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

solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
t = np.array(solfoam.TimestepValues)
N=t.size

# create a new 'Integrate Variables'
int1 = pv.IntegrateVariables(Input=solfoam)

outfile=open("presavg.dat","w")
for i in range(N):
    pv.UpdatePipeline(time=t[i], proxy=int1)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int1) )
    vol     = idat.CellData['Volume'].item()
    presint = idat.PointData['p_rgh'].item()
    print("processing time = %e\t%e\t%e\t%e" % (t[i],vol,presint,presint/vol))
    outfile.write("%e\t%e\n"%(t[i],presint/vol))

outfile.close()
