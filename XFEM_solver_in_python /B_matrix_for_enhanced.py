import numpy as np
from Determinant_function import my_funcdet
from Determinant_function_partitioning import my_funcdetforpartition

def my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,Dvalues1,
                       xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc,phivalforpart):          #To build the B-matrix for enriched elements
##phivalues added and Dvalues1 removed        
         cc=i;dd=j
         #dNdxi,detjenr,Jac = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc); 
         detjenr = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
         dNdxi,detjenrpart,Jac = my_funcdetforpartition(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,parlist,part,k);
         
       
         invj=np.linalg.inv(Jac)                #To find inverse of the jacobi matrix
         #print (detj)
        # phivalues1=phivalues[i]       ##Dummy check
         phivalues2=phivalforpart[part]     
         phivalues1=phivalues2[k]
         ##New addition(01-04-23)
         if phivalues1[j]==1:
             D=Dstd
             DD="std"
         elif phivalues1[j]==0:
             if i in elementcuthole:
                    D=D1;
                    DD="Dhole"
             else:
                    D=D2;
                    DD="Dinc"
         else:
             if i in elementcuthole:
                    D=D1;
                    DD="Dhole"
             else:
                    D=D2;
                    DD="Dinc"
             
         Benr1=np.matmul(invj,dNdxi);
         Dvalues1.append(DD);
         
        
          
          
         Benr=np.array([[(Benr1[0,0]),0,Benr1[0,0]*float(phivalues1[0]),0,Benr1[0,1],0,Benr1[0,1]*float(phivalues1[1]),0,Benr1[0,2],0,\
                      Benr1[0,2]*float(phivalues1[2]),0,Benr1[0,3],0,Benr1[0,3]*float(phivalues1[3]),0],\
                      [0,Benr1[1,0],0,Benr1[1,0]*float(phivalues1[0]),0,Benr1[1,1],0,Benr1[1,1]*float(phivalues1[1]),0,Benr1[1,2],0,Benr1[1,2]*float(phivalues1[2]),\
                       0,Benr1[1,3],0,Benr1[1,3]*float(phivalues1[3])],\
                         [Benr1[1,0],Benr1[0,0],Benr1[1,0]*float(phivalues1[0]),Benr1[0,0]*float(phivalues1[0]),Benr1[1,1],Benr1[0,1],\
                          Benr1[1,1]*float(phivalues1[1]),Benr1[0,1]*float(phivalues1[1]),Benr1[1,2],Benr1[0,2],Benr1[1,2]*float(phivalues1[2]),Benr1[0,2]*float(phivalues1[2]),\
                          Benr1[1,3],Benr1[0,3],Benr1[1,3]*float(phivalues1[3]),Benr1[0,3]*float(phivalues1[3])]]);
         return Benr,detjenr,detjenrpart,D,Dvalues1; 
        #return Benr,detjenr,D
