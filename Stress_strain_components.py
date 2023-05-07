import numpy as np

def stress_strain_component(i,disp,flatten,element,elementc,elementnoc,elementinsidehole,Dstd,D1,phivalues):
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
     if i in elementnoc:
       stress2=np.matmul(Dstd,strain2);
     elif i in elementinsidehole:
       stress2=np.matmul(D1,strain2);
  elif i in elementc:
       Benr1=np.round(np.matmul(invj,dNdxi1),2);
       Benr=np.array([[(Benr1[0,0]),0,Benr1[0,0]*float(phivalues1[0]),0,Benr1[0,1],0,Benr1[0,1]*float(phivalues1[1]),0,Benr1[0,2],0,\
                    Benr1[0,2]*float(phivalues1[2]),0,Benr1[0,3],0,Benr1[0,3]*float(phivalues1[3]),0],\
                    [0,Benr1[1,0],0,Benr1[1,0]*float(phivalues1[0]),0,Benr1[1,1],0,Benr1[1,1]*float(phivalues1[1]),0,Benr1[1,2],0,Benr1[1,2]*float(phivalues1[2]),\
                     0,Benr1[1,3],0,Benr1[1,3]*float(phivalues1[3])],\
                       [Benr1[1,0],Benr1[0,0],Benr1[1,0]*float(phivalues1[0]),Benr1[0,0]*float(phivalues1[0]),Benr1[1,1],Benr1[0,1],\
                        Benr1[1,1]*float(phivalues1[1]),Benr1[0,1]*float(phivalues1[1]),Benr1[1,2],Benr1[0,2],Benr1[1,2]*float(phivalues1[2]),Benr1[0,2]*float(phivalues1[2]),\
                        Benr1[1,3],Benr1[0,3],Benr1[1,3]*float(phivalues1[3]),Benr1[0,3]*float(phivalues1[3])]]);
       strain2=np.matmul(Benr,noddisp)
       stress2=np.matmul(D1,strain2);         ###Please change it today for better accuracy
       
       
  return stress2,strain2