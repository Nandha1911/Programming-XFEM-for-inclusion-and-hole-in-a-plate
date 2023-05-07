import numpy as np
from Determinant_function import my_funcdet

def my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc):            #To build the B-matrix for standard elements
        
         cc=i;dd=j
         dNdxi,detj,Jac = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
         #dNdxi,detj,Jac = my_funcdet(cc,dd);
         #detj=np.linalg.det(Jac)                 #To find determinant of the jacobi matrix
         invj=np.linalg.inv(Jac)                 #To find inverse of the jacobi matrix
         #c2=xi2[i];
         Bstd1=np.round(np.matmul(invj,dNdxi),2);
         
         Bstd =np.array([[Bstd1[0,0],0,Bstd1[0,1],0,Bstd1[0,2],0,Bstd1[0,3],0],\
                        [0,Bstd1[1,0],0,Bstd1[1,1],0,Bstd1[1,2],0,Bstd1[1,3]],\
                           [Bstd1[1,0],Bstd1[0,0],Bstd1[1,1],Bstd1[0,1],\
                             Bstd1[1,2],Bstd1[0,2],Bstd1[1,3],Bstd1[0,3]]]);
             
         #Bstd=np.matmul(invj,Bstd1);
         
         
         return Bstd,detj