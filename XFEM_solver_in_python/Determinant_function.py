import numpy as np

def my_funcdet(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc):   #Separate function to use determinant for both Stiffness matrix and force vector
    #dN1dxi1=N1.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);   #Differentation w.r.t c1
   # dN2dxi1=N2.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dN3dxi1=N3.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dN4dxi1=N4.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dN1dxi2=N1.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);   #Differentation w.r.t c2
    #dN2dxi2=N2.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dN3dxi2=N3.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dN4dxi2=N4.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
    #dNdxi1= np.array([[float(dN1dxi1),float(dN2dxi1),float(dN3dxi1),float(dN4dxi1)],\
    #                 [float(dN1dxi2),float(dN2dxi2),float(dN3dxi2),float(dN4dxi2)]]);
   # a=xi1[j]   #Checking
   # b=xi2[j]    Checking
    dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])  #Use this method(correct)
    dNdxi1=np.round(dNdxi1,2)
    kk=element[i]
    Jac1=np.float32(np.matmul(dNdxi1,kk));           #jacobi matrix
    detj1=np.linalg.det(Jac1);                               #To find determinant of the jacobi matrix
    
    if i in elementnoc:
         return dNdxi1,detj1,Jac1;
    else:
         return detj1
    
   # return dNdxi1,detj1,Jac1;     #Temporary
