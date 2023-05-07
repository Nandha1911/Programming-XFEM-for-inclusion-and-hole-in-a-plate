import numpy as np


def fixednodes(a,b,c,globnodes,globnumber,l):
    bdnodes1=[]
    for i in range(a):
              if c == "Edge1" :
                  j=[round(i*l,2),0];
                  k=globnodes.index(j);
                  l1=np.array([k,k+1]);
              elif c == "Edge2" :
                  j=[round(b*l,2),round(i*l,2)];
                  k=globnodes.index(j);
                  l1=np.array([k+(b+1),k]);
              elif c == "Edge3":
                  j=[round(i*l,2),round(b*l,2)]; 
                  k=globnodes.index(j);
                  l1=np.array([k,k+1]);
              elif c == "Edge4":
                  j=[0,round(i*l,2)];
                  k=globnodes.index(j);
                  l1=np.array([k,k+(b+1)]);
              #for z in range(len(globnumber)):
                  #if ((globnumber[z]==l1).all()):
              bdnodes1.append(l1);
              bdnodes2=np.unique(bdnodes1)
    bdnodes=bdnodes2.tolist()
    return bdnodes