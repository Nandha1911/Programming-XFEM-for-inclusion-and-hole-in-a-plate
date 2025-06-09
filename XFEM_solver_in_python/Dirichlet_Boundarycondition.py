import numpy as np

def DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3,F):
  # for j in range(len(fixedDOFSlist[i])):
      l1=np.array(fixDOF)
      l2=flatten1.index(l1)
      l3.append(l2)
      for m in range(len(flatten1)):   ###For column side
          if m == l2:
              K[m,l2]=1
          else:
              K[m,l2]=0
              K[l2,m]=0
      for n in range(len(flatten1)):   ###For row side
          if n == l2:
               F[l2,0]=0
          #else:
             #  K[l2,n]=0
    
      return K,l3,F
