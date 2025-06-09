import numpy as np


def Partitioning_method(i,l,elementc,element):
  norforpartition=3
  nocforpartition=3
  element2=[];
  le=float(format((l/norforpartition),".3f"));
  for m in range(norforpartition):                  ####Represents columns (3*3=Total 9 elements)
    le1=float(format((l/norforpartition),".3f"));
    k=elementc[i]
    p=element[k]
    for n in range(nocforpartition):                 ####Represents rows   (3*3=Total 9 elements)
              
             # node3=([p[1][1]+(le1),p[1][0]+(le*m)]);
             # node2=([p[0][0]+(le1),p[0][1]+(le*m)]);
             # node2=([p[1][0]+(le*n),p[1][1]+(le*m)]);
              #node2=([p[0][0]+(le*n),p[0][1]+(le*(m+1))]);
              node1=([p[1][0]+(le*(n)),p[1][1]+(le*(m+1))]);
              node2=([p[1][0]+(le*n),p[1][1]+(le*m)]);
              node3=([p[1][0]+(le1),p[1][1]+(le*m)]);
              node4=([p[1][0]+(le1),p[1][1]+(le*(m+1))]);
              #node4=([p[0][0]+(le1),p[0][1]+(le*m)]);
              #node4=([p[1][0]+(le*n),p[1][1]+(le*(m+1))]);
              #if m==0 and n==2:
                #  node4=(p[3][0],p[3][1])
              #elif m==2 and n==0:
                 # node2=(p[1][0],p[1][1])
              #elif m==2 and n==2:
                #  node3=(p[2][0],p[2][1])
                  
              element1=np.array([node1,node2,node3,node4]);
              for t in range(len(element1)):
                  s=element1[t];
                  if (p[3][0]-s[0]) <= 0.01 and (p[3][1]-s[1]) <= 0.01:
                      s[0]=p[3][0]
                      s[1]=p[3][1]
                      element1[t]=(s[0],s[1])
                  elif (p[3][0]-s[0]) <= 0.01:
                      s[0]=p[3][0]
                      element1[t]=(s[0],s[1])
                  elif (p[3][1]-s[1]) <= 0.01:
                      s[1]=p[3][1]
                      element1[t]=(s[0],s[1])
                  #elif (p[2][0]-s[0]) <= 0.02 and (p[2][1]-s[1]) <= 0.02:
                      #s[0]=p[2][0]
                      #s[1]=p[2][1]
                     # element1[t]=(s[0],s[1])
                  #element1=np.array([node1,node2,node3,node4]);
              element2.append(element1)
              le1=le1+le
              
  return element2
