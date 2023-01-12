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


lf=2.25 #mm
hc=0.5  #mm
spacrad=0.125 #mm
Db=1.6e-9
vp=float(argv[3])

dh=4.0*(lf*hc-np.pi*spacrad**2)/(2.0*lf+2*np.pi*spacrad)
print(dh)

dhbyD=dh*1e-3/Db

solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['membrane1']
solfoam.CellArrays = ['C']
cb=float(argv[1])
rs=float(argv[2])
t = np.array(solfoam.TimestepValues)
N=t.size

# Properties modified on calculator1
calculator1 = pv.Calculator(Input=solfoam)
calculator1.AttributeType = 'Point Data'
calculator1.ResultArrayName = 'cp'
calculator1.Function = 'C/%5.2e'%(cb)

calculator2 = pv.Calculator(Input=calculator1)
calculator2.AttributeType = 'Point Data'
calculator2.ResultArrayName = 'k'
calculator2.Function = '%5.2e*%5.2e/ln(%5.2e/(1.0/cp+%5.2e-1.0))'%(dhbyD,vp,rs,rs)
#calculator2.Function = '%5.2e*%5.2e'%(dhbyD,vp)
print(calculator2.Function)

# create a new 'Integrate Variables'
int2 = pv.IntegrateVariables(Input=calculator2)

outfile=open("cpavg.dat","w")
for i in range(N):
    pv.UpdatePipeline(time=t[i], proxy=int2)
    idat    = dsa.WrapDataObject(pv.servermanager.Fetch(int2) )
    area    = idat.CellData['Area'].item()
    cpint   =idat.PointData['cp'].item()
    kint    =idat.PointData['k'].item()
    Cint    =idat.PointData['C'].item()
    print("processing time = %e\t%e\t%e\t%e\t%e\t%e" % (t[i],area,cpint,cpint/area,kint/area,Cint/area))
    outfile.write("%e\t%e\t%e\t%e\n"%(t[i],cpint/area,kint/area,Cint/area))

outfile.close()
