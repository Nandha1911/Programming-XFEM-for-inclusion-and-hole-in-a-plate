import numpy as np


def Tractionforce(a,b,c,globnodes,globnumber,l):
    bdelement=[]
    check=[]
    for i in range(a+1):
              if c == "Edge1" :
                  j=np.array([round(i*l,2),0]);
                  for k in range(len(globnodes)):
                      check1=np.allclose(j,globnodes[k])
                      check.append(check1);
                      if check1==True:
                          break;
                  l1=np.array([k+(a+1),k,k+1,k+(a+2)]);
              elif c == "Edge2" :
                  j=np.array([round(b*l,2),round(i*l,2)]);
                  for k in range(len(globnodes)):
                      check1=np.allclose(j,globnodes[k])
                      check.append(check1);
                      if check1==True:
                          break;
                  l1=np.array([k+(b),k-1,k,k+(b+1)]);
              elif c == "Edge3":
                  j=np.array([round(i*l,2),round(b*l,2)]); 
                  for k in range(len(globnodes)):
                      check1=np.allclose(j,globnodes[k])
                      check.append(check1);
                      if check1==True:
                          break;
                  l1=np.array([k,k-1-(a),k-(a),k+1]);
              elif c == "Edge4":
                  j=np.array([0,round(i*l,2)]); 
                  for k in range(len(globnodes)):
                      check1=np.allclose(j,globnodes[k])
                      check.append(check1);
                      if check1==True:
                          break;
                  l1=np.array([k+(b+1),k,k+1,k+(b+2)]);
              for z in range(len(globnumber)):
                  if ((globnumber[z]==l1).all()):
                      bdelement.append(z);
                      break;
    return bdelement
