import numpy as np

def creatingelement(i,j,l):

     node1=([(j*l),(i*l)]);
     node2=([(j*l),(i+1)*l]);
     node3=([(j+1)*l,(i+1)*l]);
     node4=([(j+1)*l,i*l]);

     element1=np.array([node1,node2,node3,node4]);
     return element1;