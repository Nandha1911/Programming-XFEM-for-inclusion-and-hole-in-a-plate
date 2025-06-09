import numpy as np

def creatingelement(i,j,l,noc,nor):

     node1=([(j*l),(i+1)*l]);
     node2=([(j*l),(i*l)]);
     node3=([(j+1)*l,i*l]);
     node4=([(j+1)*l,(i+1)*l]);

     element1=np.array([node1,node2,node3,node4]);
     if i ==0 and j==0:
         cornerelements1=1;
     elif i==0 and j==nor-1:
         cornerelements1=1;
     elif i==noc-1 and j==0:
         cornerelements1=1;
     elif i==noc-1 and j==nor-1:
         cornerelements1=1;
     else:
         cornerelements1=0;
     return element1,cornerelements1;
