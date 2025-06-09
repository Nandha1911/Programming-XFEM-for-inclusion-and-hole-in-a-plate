import numpy as np

def Tottractforceassm(i,globnumber,elementnoc,DOF,globnodesextra,fbloc2,Ftract1):
        #Ftract1=np.zeros(np.shape(flatten));
        #fbloc1=totbodyforce[i];
        if i in elementnoc:
            k=globnumber[i];
            mapcol=[];
            for j in range(len(globnumber[i])): #To calculate for Standard elements
                map2=[];
                for d in range(DOF):
                   map1=int((((k[j]-1)*DOF)+(d+1))+1+np.sum(globnodesextra[0:k[j]]));
                   map2.append(map1);
                mapcol.append(map2);
            #K=np.zeros((totDOF,totDOF));
           # print(len(map2));
          #  for m in range(len(map2)-1):
            for n in range(len(mapcol)):
                  fbglob1=fbloc2[(n*2):((n*2)+2)];
                  fbglob2=Ftract1[mapcol[n][0]:(mapcol[n][1]+1)];
                  fbglob=np.add(fbglob1,fbglob2);
                  Ftract1[mapcol[n][0]:(mapcol[n][1]+1)]=fbglob;
          
          
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
                     # for m in range(len(map2)-1):
                      for n in range(len(mapcol)):
                            # Kglob=Kloc1[0:2,0:2]+K[mapcol[0,0]:(mapcol[0,1]+1)];
                            #K=np.zeros((totDOF,totDOF));
                            fbglob1=fbloc2[(n*4):((n*4)+4)];
                            fbglob2=Ftract1[mapcol[n][0]:(mapcol[n][3]+1)];
                            fbglob=np.add(fbglob1,fbglob2);
                            Ftract1[mapcol[n][0]:(mapcol[n][3]+1)]=fbglob;
        
        return Ftract1
