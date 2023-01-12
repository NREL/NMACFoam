import numpy as np
from feed_channel_data import *

def write_ofoam_preamble(outfile):

    outfile.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
    outfile.write("| =========                 |                                                 |\n")
    outfile.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
    outfile.write("|  \\    /   O peration     | Version:  5                                     |\n")
    outfile.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
    outfile.write("|    \\/     M anipulation  |                                                 |\n")
    outfile.write("\*---------------------------------------------------------------------------*/\n")
    outfile.write("FoamFile\n")
    outfile.write("{\n")
    outfile.write("\tversion     2.0;\n")
    outfile.write("\tformat      ascii;\n")
    outfile.write("\tclass       dictionary;\n")
    outfile.write("\tobject      blockMeshDict;\n")
    outfile.write("}\n\n")
    outfile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")
    outfile.write("convertToMeters 0.001;\n\n")
    outfile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

def write_vertices(outfile):

    outfile.write("\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    outfile.write("vertices\n(\n")

    counter=0

    for yiter in range(2):
        for ziter in range(2):
    
            outfile.write("(0.0 "+str(yiter*channel_span)+" "+\
                    str(ziter*channel_W)+ ") // "+str(counter)+"\n")
            counter=counter+1

            offset=sp_off1;
            for sp in range(nspacers):
                spac_center=offset+sp_gap*sp;
                outfile.write("( "+str(spac_center-0.5*channel_W)+" "+\
                        str(yiter*channel_span)+" "+str(ziter*channel_W) + ") // "+str(counter)+"\n")
                counter=counter+1
                outfile.write("( "+str(spac_center+0.5*channel_W)+" "+\
                        str(yiter*channel_span)+" "+str(ziter*channel_W) + ") // "+str(counter)+"\n")
                counter=counter+1
            
            outfile.write("("+str(channel_L)+"  "+str(yiter*channel_span)+" "+\
                    str(ziter*channel_W)+ ") // "+str(counter)+"\n")
            counter=counter+1

            outfile.write("\n\n");

    #circle points
    for yiter in range(2):
        for ziter in range(2):
            
            ziter_mod=2*ziter-1 #in -1,1 space
            offset=sp_off1;
            for sp in range(nspacers):
                spac_center=offset+sp_gap*sp;
                outfile.write("( "+str(spac_center-sp_R_x)+" "+\
                        str(yiter*channel_span)+" "+str(0.5*channel_W+ziter_mod*sp_R_z) + ") // "+str(counter)+"\n")
                counter=counter+1
                outfile.write("( "+str(spac_center+sp_R_x)+" "+\
                        str(yiter*channel_span)+" "+str(0.5*channel_W+ziter_mod*sp_R_z) + ") // "+str(counter)+"\n")
                counter=counter+1
            
            outfile.write("\n\n");

    outfile.write(");\n")

def write_edges(outfile):

    outfile.write("\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    outfile.write("edges\n(\n")
    
    ncart_points=(2*nspacers+2)*2*2

    offset=sp_off1;
    for yiter in range(2):
        for sp in range(nspacers):

            spac_center=offset+sp_gap*sp;
            ind1=yiter*nspacers*4+ncart_points+sp*2
            ind2=ind1+1
            ind3=yiter*nspacers*4+ncart_points+nspacers*2+sp*2
            ind4=ind3+1

            outfile.write("arc "+str(ind1)+" "+str(ind2)+" ")
            outfile.write("( "+str(spac_center)+" "+str(yiter*channel_span)+" "+str(0.5*channel_W-sp_R)+" )\n")  

            outfile.write("arc "+str(ind3)+" "+str(ind4)+" ")
            outfile.write("( "+str(spac_center)+" "+str(yiter*channel_span)+" "+str(0.5*channel_W+sp_R)+" )\n")  
        
            outfile.write("arc "+str(ind1)+" "+str(ind3)+" ")
            outfile.write("( "+str(spac_center-sp_R)+" "+str(yiter*channel_span)+" "+str(0.5*channel_W)+" )\n")  
        
            outfile.write("arc "+str(ind2)+" "+str(ind4)+" ")
            outfile.write("( "+str(spac_center+sp_R)+" "+str(yiter*channel_span)+" "+str(0.5*channel_W)+" )\n")  

            outfile.write("\n\n")
        

    outfile.write(");\n")

def write_this_block(outfile,comment,ids,mesh,zonename="none"):
    
    outfile.write("\n //"+comment+"\n")
    outfile.write("hex (")
    for i in range(len(ids)):
        outfile.write(str(ids[i])+" ")
    outfile.write(")\n")

    if(zonename != "none"):
        outfile.write(zonename+"\n")

    outfile.write("( %d %d %d )\n"%(mesh[0],mesh[1],mesh[2]))
    outfile.write("SimpleGrading (1 1 1)\n")


def write_blocks(outfile):

    outfile.write("\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    outfile.write("blocks\n(\n")

    idarray   = np.zeros(8,dtype=np.int)
    mesharray = np.zeros(3,dtype=np.int)

    npts_along_line=2*nspacers+2
    nxz_pts=npts_along_line*2

    blocknum=0
    #cartesian blocks
    for cblocks in range(int(npts_along_line/2)):
        
        if(cblocks==0):
            blocklen=sp_off1-0.5*channel_W
        elif(cblocks==nspacers):
            blocklen=channel_L-(sp_off1+(nspacers-1)*sp_gap)
        else:
            blocklen=sp_gap-channel_W

        idarray[0] = 2*cblocks
        idarray[1] = npts_along_line+2*cblocks
        idarray[2] = npts_along_line+2*cblocks+1
        idarray[3] = 2*cblocks+1

        idarray[4] = nxz_pts+idarray[0];
        idarray[5] = nxz_pts+idarray[1];
        idarray[6] = nxz_pts+idarray[2];
        idarray[7] = nxz_pts+idarray[3];
        
        mesharray[0] = int(nz*channel_W)+2
        mesharray[1] = int(nx*blocklen)+2
        mesharray[2] = ny
        zonename="none"
        
        write_this_block(outfile,"block %d"%(blocknum),idarray,mesharray,zonename)
        blocknum=blocknum+1


    ncart_points=(2*nspacers+2)*2*2

    #spacer blocks
    for sp in range(nspacers):
        
        idarray[0] = 2*sp+1
        idarray[1] = npts_along_line+2*sp+1
        idarray[2] = ncart_points+2*nspacers+2*sp
        idarray[3] = ncart_points+2*sp

        idarray[4] = nxz_pts+idarray[0];
        idarray[5] = nxz_pts+idarray[1];
        idarray[6] = 4*nspacers+idarray[2];
        idarray[7] = 4*nspacers+idarray[3];
        
        mesharray[0] = int(nz*channel_W)+2
        mesharray[1] = int(nr*(channel_W/np.sqrt(2)-sp_R))+2
        mesharray[2] = ny
        zonename="none"
        
        write_this_block(outfile,"block %d"%(blocknum),idarray,mesharray,zonename)
        blocknum=blocknum+1
        
        idarray[0] = 2*sp+1
        idarray[1] = ncart_points+2*sp
        idarray[2] = idarray[1]+1
        idarray[3] = idarray[0]+1

        idarray[4] = nxz_pts+idarray[0];
        idarray[5] = 4*nspacers+idarray[1];
        idarray[6] = 4*nspacers+idarray[2];
        idarray[7] = nxz_pts+idarray[3];
        
        mesharray[0] = int(nr*(channel_W/np.sqrt(2)-sp_R))+2
        mesharray[1] = int(nx*channel_W)+2
        mesharray[2] = ny
        zonename="none"
        
        write_this_block(outfile,"block %d"%(blocknum),idarray,mesharray,zonename)
        blocknum=blocknum+1
        
        idarray[0] = ncart_points+2*nspacers+2*sp+1
        idarray[1] = npts_along_line+2*sp+2
        idarray[2] = idarray[1]-npts_along_line
        idarray[3] = idarray[0]-2*nspacers

        idarray[4] = 4*nspacers+idarray[0];
        idarray[5] = nxz_pts+idarray[1];
        idarray[6] = nxz_pts+idarray[2];
        idarray[7] = 4*nspacers+idarray[3];
        
        mesharray[0] = int(nr*(channel_W/np.sqrt(2)-sp_R))+2
        mesharray[1] = int(nz*channel_W)+2
        mesharray[2] = ny
        zonename="none"
        
        write_this_block(outfile,"block %d"%(blocknum),idarray,mesharray,zonename)
        blocknum=blocknum+1
        
        idarray[0] = ncart_points+2*nspacers+2*sp
        idarray[1] = npts_along_line+2*sp+1
        idarray[2] = idarray[1]+1
        idarray[3] = idarray[0]+1

        idarray[4] = 4*nspacers+idarray[0];
        idarray[5] = nxz_pts+idarray[1];
        idarray[6] = nxz_pts+idarray[2];
        idarray[7] = 4*nspacers+idarray[3];
        
        mesharray[0] = int(nr*(channel_W/np.sqrt(2)-sp_R))+2
        mesharray[1] = int(nx*channel_W)+2
        mesharray[2] = ny
        zonename="none"
        
        write_this_block(outfile,"block %d"%(blocknum),idarray,mesharray,zonename)
        blocknum=blocknum+1
        
    outfile.write(");\n")




def write_patches(outfile):

    outfile.write("\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
    outfile.write("boundary\n(\n")
    
    npts_along_line=2*nspacers+2
    nxz_pts=npts_along_line*2
    ncart_points=(2*nspacers+2)*2*2
    
    outfile.write("\n\tinlet\n\t{\n")

    outfile.write("\t\ttype patch;\n\t\tfaces\n\t\t(\n")
    ind1=0
    ind2=ind1+nxz_pts;
    ind3=ind2+npts_along_line;
    ind4=ind1+npts_along_line;

    outfile.write("\t\t( ")
    outfile.write(str(ind1)+" ")
    outfile.write(str(ind2)+" ")
    outfile.write(str(ind3)+" ")
    outfile.write(str(ind4)+")\n")
    
    outfile.write("\t\t); \n")
    outfile.write("\t} \n")

    outfile.write("\n\toutlet\n\t{\n")
    
    ind1=npts_along_line-1
    ind2=ind1+nxz_pts
    ind3=ind2+npts_along_line
    ind4=ind1+npts_along_line

    outfile.write("\t\ttype patch;\n\t\tfaces\n\t\t(\n")
    outfile.write("\t\t( ")
    outfile.write(str(ind1)+" ")
    outfile.write(str(ind2)+" ")
    outfile.write(str(ind3)+" ")
    outfile.write(str(ind4)+")\n")
    
    outfile.write("\t\t); \n")
    outfile.write("\t} \n")


    outfile.write("\n\tmembrane1\n\t{\n")
    
    outfile.write("\t\ttype patch;\n\t\tfaces\n\t\t(\n")

    for i in range(npts_along_line-1):

        ind1=i
        ind2=ind1+nxz_pts
        ind3=ind2+1
        ind4=ind1+1

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")
    
    outfile.write("\t\t); \n")
    outfile.write("\t} \n")
    
    outfile.write("\n\tmembrane2\n\t{\n")
    outfile.write("\t\ttype patch;\n\t\tfaces\n\t\t(\n")

    for i in range(npts_along_line-1):

        ind1=i+npts_along_line
        ind2=ind1+nxz_pts
        ind3=ind2+1
        ind4=ind1+1

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")
    
    outfile.write("\t\t); \n")
    outfile.write("\t} \n")
    
    outfile.write("\n\tspacers\n\t{\n")
    outfile.write("\t\ttype patch;\n\t\tfaces\n\t\t(\n")

    for sp in range(nspacers):

        ind1=ncart_points+2*sp
        ind2=ind1+4*nspacers
        ind3=ind2+1
        ind4=ind1+1

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")
        
        ind1=ncart_points+2*nspacers+2*sp
        ind2=ind1+4*nspacers
        ind3=ind2+1
        ind4=ind1+1

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")
        
        ind1=ncart_points+2*sp
        ind2=ind1+2*nspacers
        ind3=ind2+4*nspacers
        ind4=ind1+4*nspacers

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")
        
        ind1=ncart_points+2*sp+1
        ind2=ind1+2*nspacers
        ind3=ind2+4*nspacers
        ind4=ind1+4*nspacers

        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        outfile.write("\n\n")
    
    outfile.write("\t\t); \n")
    outfile.write("\t} \n")
    
    outfile.write("\n\tfront\n\t{\n")
    outfile.write("\t\ttype empty;\n\t\tfaces\n\t\t(\n")
    
    #cartesian blocks
    for cblocks in range(int(npts_along_line/2)):
        
        ind1 = 2*cblocks
        ind2 = npts_along_line+2*cblocks
        ind3 = npts_along_line+2*cblocks+1
        ind4 = 2*cblocks+1
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

    #spacer blocks
    for sp in range(nspacers):
        
        ind1 = 2*sp+1
        ind2 = npts_along_line+2*sp+1
        ind3 = ncart_points+2*nspacers+2*sp
        ind4 = ncart_points+2*sp
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = 2*sp+1
        ind2 = ncart_points+2*sp
        ind3 = ind2+1
        ind4 = ind1+1
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = ncart_points+2*nspacers+2*sp+1
        ind2 = npts_along_line+2*sp+2
        ind3 = ind2-npts_along_line
        ind4 = ind1-2*nspacers
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = ncart_points+2*nspacers+2*sp
        ind2 = npts_along_line+2*sp+1
        ind3 = ind2+1
        ind4 = ind1+1
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

    outfile.write("\t\t); \n")
    outfile.write("\t} \n")
    
    outfile.write("\n\tback\n\t{\n")
    outfile.write("\t\ttype empty;\n\t\tfaces\n\t\t(\n")
    
    #cartesian blocks
    for cblocks in range(int(npts_along_line/2)):
        
        ind1 = 2*cblocks+nxz_pts
        ind2 = npts_along_line+2*cblocks+nxz_pts
        ind3 = npts_along_line+2*cblocks+1+nxz_pts
        ind4 = 2*cblocks+1+nxz_pts
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

    #spacer blocks
    for sp in range(nspacers):
        
        ind1 = 2*sp+1+nxz_pts
        ind2 = npts_along_line+2*sp+1+nxz_pts
        ind3 = ncart_points+2*nspacers+2*sp+4*nspacers
        ind4 = ncart_points+2*sp+4*nspacers
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = 2*sp+1+nxz_pts
        ind2 = ncart_points+2*sp+4*nspacers
        ind3 = ind2+1
        ind4 = ind1+1
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = ncart_points+2*nspacers+2*sp+1+4*nspacers
        ind2 = npts_along_line+2*sp+2+nxz_pts
        ind3 = ind2-npts_along_line
        ind4 = ind1-2*nspacers
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

        ind1 = ncart_points+2*nspacers+2*sp+4*nspacers
        ind2 = npts_along_line+2*sp+1+nxz_pts
        ind3 = ind2+1
        ind4 = ind1+1
        
        outfile.write("\t\t( ")
        outfile.write(str(ind1)+" ")
        outfile.write(str(ind2)+" ")
        outfile.write(str(ind3)+" ")
        outfile.write(str(ind4)+")\n")

    outfile.write("\t\t); \n")
    outfile.write("\t} \n")


    outfile.write(");\n")

    

#main
outfile=open("blockMeshDict_channel","w")
write_ofoam_preamble(outfile)
write_vertices(outfile)
write_edges(outfile)
write_blocks(outfile)
write_patches(outfile)
outfile.close()
