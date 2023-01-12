# NMACFoam
Numerical-model for Membrane-spacer Assemblies and Configurations in OpenFOAM

# Build instructions
* The solver is now compatible with openfoam version 6 - see https://github.com/OpenFOAM/OpenFOAM-dev/tree/version-6 
and will be upgraded to the latest openfoam soon.
* source openfoam's bashrc - source <path-to-OpenFoam-dev>/etc/bashrc
* go to solver folder (cd NMACFoam) and do $wmake
  
# Test cases
  Go to each of the following test cases and use the Allrun script which essentially creates a mesh and runs the solver
  You can run in parallel as well after decomposition using decomposePar and using mpirun -n nprocs NMACFoam -parallel
  
* ## flatMembrane
  flow through channel with one boundary being a membrane (follows the paper by A. Subramani et al. Journal of Membrane Science 277 (2006) 7–17)
  ![flatMem](https://user-images.githubusercontent.com/7399475/119549154-2e806100-bd54-11eb-94cd-3667268006ee.png)
 
* ## membraneSpacer
  flow through channel with spacers (follows the paper by Guillan and Hoek, Chemical Engineering Journal 149 (2009) 221–231)
 
  ![memspacer](https://user-images.githubusercontent.com/7399475/119550594-d64a5e80-bd55-11eb-8e84-33c858c88e77.png)
  ![NMAC_GuilanHoek_pressure_compare](https://user-images.githubusercontent.com/7399475/120541953-c3fca000-c3a7-11eb-9fe5-9ab4a88cc60f.png)
  ![NMAC_GuilanHoek_Shcorr](https://user-images.githubusercontent.com/7399475/120541968-c6f79080-c3a7-11eb-9b2b-aa6287d348d6.png)

 * ## one/twoPorousWalls
  verification case for pressure drop for flow through channel with porous walls
  
  ![pressuredrop](https://user-images.githubusercontent.com/7399475/119553401-e0ba2780-bd58-11eb-9c35-b82ac01c2a50.png)

  
