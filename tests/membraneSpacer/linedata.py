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
#solfoam.CellArrays = ['C']
solfoam.PointArrays = ['C','p']
t = np.array(solfoam.TimestepValues)
N=t.size

# create a new 'Cell Data to Point Data'
#cellDatatoPointData1 = pv.CellDatatoPointData(Input=solfoam)
#cellDatatoPointData1.CellDataArraytoprocess = ['C']

# create a new 'Slice'
#slice1 = pv.Slice(Input=cellDatatoPointData1)

pltline = pv.PlotOverLine(Input=solfoam,
            Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
pltline.Source.Point1 = [0.0,   0.0025, 1e-8]
pltline.Source.Point2 = [0.035, 0.0025, 1e-8]


pv.UpdatePipeline(time=t[-1], proxy=pltline)
idat    = dsa.WrapDataObject(pv.servermanager.Fetch(pltline))
outarr=np.array([idat.Points[:,0],idat.PointData['C'],idat.PointData['p']])
np.savetxt("conc_and_p.txt", np.transpose(outarr), delimiter=" ")


