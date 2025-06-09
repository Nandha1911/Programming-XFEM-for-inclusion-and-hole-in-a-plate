import numpy as np

def Disp_stress_strain(i,disp,flatten,sigma,Bmatrix,elementc,elementnoc,elementinsidehole,phivalues):
    disp3=[];
    for j in range(len(disp)):
        disp1=disp[j]
        for k in range(len(disp1)):
            #call=disp1[k]
            disp2=flatten[disp1[k]]
            disp3.append(disp2)
    noddisp=np.array(disp3);          ###Nodal displacements of a single element
    c1 = 0;                #Since Displacement interpolation is taken at element centroid,it is taken as zero
    c2 = 0;
    N1=((1-c1)*(1-c2))/4;
    N2=((1+c1)*(1-c2))/4;      
    N3=((1+c1)*(1+c2))/4;
    N4=((1-c1)*(1+c2))/4;
    c3=np.array([[N1],[N2],[N3],[N4]])
    if i in elementnoc:
         cx=np.array([noddisp[0],noddisp[2],noddisp[4],noddisp[6]])
         cy=np.array([noddisp[1],noddisp[3],noddisp[5],noddisp[7]])
         dispx=np.matmul(np.transpose(c3),(cx));
         dispy=np.matmul(np.transpose(c3),(cy));
         actdisp=[dispx,dispy]
    elif i in elementc:
         phi=phivalues[i]
         cx=np.array([noddisp[0],noddisp[4],noddisp[8],noddisp[12]])
         cy=np.array([noddisp[1],noddisp[5],noddisp[9],noddisp[13]])
         ax=np.array([noddisp[2],noddisp[6],noddisp[10],noddisp[14]])
         ay=np.array([noddisp[3],noddisp[7],noddisp[11],noddisp[15]])
         Fphi=abs(N1*phi[0]+N2*phi[1]+N3*phi[2]+N4*phi[3]);
         dispx=np.matmul(np.transpose(c3),(cx))+(Fphi*np.matmul(np.transpose(c3),(ax)));
         dispy=np.matmul(np.transpose(c3),(cy))+(Fphi*np.matmul(np.transpose(c3),(ay)));
         actdisp=[dispx,dispy]
    #displacement.append(disp3);
    #stresss1=np.matmul((sigma[i]),disp4);
    #strains1=np.matmul(Bmatrix[i],disp4);
    return disp3,actdisp#stresss1,strains1
