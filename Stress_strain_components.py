import numpy as np ###Last update+18/05--05:08

def stress_strain_component(umatrix,flatten,element,elementc,elementnoc,elementinsidehole,Dstd,D1,phivalues,\
                            phivalforpart,increment):
 stress=[];strain=[];  
 for i in range(len(element)):
   disp=umatrix[i]    
   disp3=[];
   for j in range(len(disp)):
        disp1=disp[j]
        for k in range(len(disp1)):
            #call=disp1[k]
            disp2=flatten[disp1[k]]
            disp3.append(disp2)
   noddisp=np.array(disp3);
   phivalues1=phivalues[i];   
   xi1=0;xi2=0;
   dNdxi1=np.array([[(xi2-1)/4,(1-xi2)/4,(1+xi2)/4,(-xi2-1)/4],[(xi1-1)/4,(-xi1-1)/4,(1+xi1)/4,(1-xi1)/4]])  #Use this method(correct)
   dNdxi1=np.round(dNdxi1,2)
   kk=element[i];
   Jac1=np.float32(np.matmul(dNdxi1,kk)); 
   invj=np.linalg.inv(Jac1) 
   if i in elementnoc:
     Bstd1=np.round(np.matmul(invj,dNdxi1),2);
 
     Bstd =np.array([[Bstd1[0,0],0,Bstd1[0,1],0,Bstd1[0,2],0,Bstd1[0,3],0],\
                [0,Bstd1[1,0],0,Bstd1[1,1],0,Bstd1[1,2],0,Bstd1[1,3]],\
                   [Bstd1[1,0],Bstd1[0,0],Bstd1[1,1],Bstd1[0,1],\
                     Bstd1[1,2],Bstd1[0,2],Bstd1[1,3],Bstd1[0,3]]]);
     strain2=np.matmul(Bstd,noddisp)
     if i in elementinsidehole:
       stress2=np.matmul(D1,strain2);
     else:
       stress2=np.matmul(Dstd,strain2);
   elif i in elementc:
       #phivaluesforsubgrid=phivalforpart[increment];
       Benr1=np.matmul(invj,dNdxi1)
       Benr=np.array([[(Benr1[0,0]),0,Benr1[0,0]*float(phivalues1[0]),0,Benr1[0,1],0,Benr1[0,1]*float(phivalues1[1]),0,Benr1[0,2],0,\
                    Benr1[0,2]*float(phivalues1[2]),0,Benr1[0,3],0,Benr1[0,3]*float(phivalues1[3]),0],\
                    [0,Benr1[1,0],0,Benr1[1,0]*float(phivalues1[0]),0,Benr1[1,1],0,Benr1[1,1]*float(phivalues1[1]),0,Benr1[1,2],0,Benr1[1,2]*float(phivalues1[2]),\
                     0,Benr1[1,3],0,Benr1[1,3]*float(phivalues1[3])],\
                       [Benr1[1,0],Benr1[0,0],Benr1[1,0]*float(phivalues1[0]),Benr1[0,0]*float(phivalues1[0]),Benr1[1,1],Benr1[0,1],\
                        Benr1[1,1]*float(phivalues1[1]),Benr1[0,1]*float(phivalues1[1]),Benr1[1,2],Benr1[0,2],Benr1[1,2]*float(phivalues1[2]),Benr1[0,2]*float(phivalues1[2]),\
                        Benr1[1,3],Benr1[0,3],Benr1[1,3]*float(phivalues1[3]),Benr1[0,3]*float(phivalues1[3])]]);
       if sum(phivalues1) >= 0:    
          strain2=np.matmul(Benr,noddisp)
          stress2=np.matmul(Dstd,strain2);         
       else:
          strain2=np.matmul(Benr,noddisp)
          stress2=np.matmul(D1,strain2); 
          increment=increment+1
   stress.append(stress2)
   strain.append(strain2)
 return stress,strain,increment
