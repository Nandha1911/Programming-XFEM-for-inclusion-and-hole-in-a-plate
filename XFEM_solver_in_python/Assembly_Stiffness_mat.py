import numpy as np


def assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix):
  for i in range(len(element)):
      Kloc1=Kloc[i];
      if i in elementnoc:
          k=globnumber[i];
          mapcol=[];
          for j in range(len(globnumber[i])): #To calculate for Standard elements
              map2=[];
              for d in range(DOF):
                 map1=int((((k[j]-1)*DOF)+(d+1))+1+np.sum(globnodesextra[0:k[j]]));
                 map2.append(map1);
              mapcol.append(map2);
              #umatrix.append(mapcol)
          #K=np.zeros((totDOF,totDOF));
          for m in range(len(mapcol)):
              for n in range(len(mapcol)):
                # Kglob=Kloc1[0:2,0:2]+K[mapcol[0,0]:(mapcol[0,1]+1)];
                #K=np.zeros((totDOF,totDOF));
                Kglob1=Kloc1[(m*2):((m*2)+2),(n*2):((n*2)+2)];
                Kglob2=K[mapcol[m][0]:(mapcol[m][1]+1),mapcol[n][0]:(mapcol[n][1]+1)];
                Kglob=np.add(Kglob1,Kglob2);
                K[mapcol[m][0]:(mapcol[m][1]+1),mapcol[n][0]:(mapcol[n][1]+1)]=Kglob;
        
      else:
                    k=globnumber[i];
                    mapcol=[];
                    for j in range(len(globnumber[i])): #To calculate for enriched elements
                        map2=[];
                        for d in range(DOF+2):
                           map1=int((((k[j]-1)*DOF)+(d+1))+1+np.sum(globnodesextra[0:k[j]]));
                           map2.append(map1);
                        mapcol.append(map2);
                    #K=np.zeros((totDOF,totDOF));
                    for m in range(len(mapcol)):
                        for n in range(len(mapcol)):
                          # Kglob=Kloc1[0:2,0:2]+K[mapcol[0,0]:(mapcol[0,1]+1)];
                          #K=np.zeros((totDOF,totDOF));
                          Kglob1=Kloc1[(m*4):((m*4)+4),(n*4):((n*4)+4)];
                          Kglob2=K[mapcol[m][0]:(mapcol[m][3]+1),mapcol[n][0]:(mapcol[n][3]+1)];
                          Kglob=np.add(Kglob1,Kglob2);
                          K[mapcol[m][0]:(mapcol[m][3]+1),mapcol[n][0]:(mapcol[n][3]+1)]=Kglob;
      
              
      umatrix.append(mapcol);
      
  return K,umatrix
