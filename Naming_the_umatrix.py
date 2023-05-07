import numpy as np


def naming_umatrix(noc,nor,l,new,zz,aDOFsiz,usize,globnodesextra,globnodes,globU,u,aDOF):
  for i in range(noc+1):
      for j in range(nor+1):
        check=[];             #Creating temporary list to store Boolean values in order to
        tempU=[];                      #check whether the global node array matches with any of the 
                              #new values(List) --Comparing global node value with new(variable List
        globnode=[round(j*l,2),round(i*l,2)];        #problem in globnode--0.1cm mesh
        for k in range(len(new)):  #To check the new(varible) list,we are using for loop(k): 
            check1=np.array_equal(globnode,new[k]);
            check.append(check1);
            if check1==True:
                break;
        if any(check) == True:
                u1='x'+str(zz);
                u2='y'+str(zz);
                a1='a'+str(zz);
                a2='b'+str(zz);
                uu=np.array([u1,u2]);
                rec=np.array([a1,a2]);
                aDOFsiz=aDOFsiz+2;           #Use to Build Assembly matrix 
                usize=usize+2;               #Use to Build Assembly matrix
                u.append(uu);
                
                aDOF.append(rec);
                globnodes1=2;
                #globnodesextra.append(globnodes1);
                
            
        else:    
               u1='x'+str(zz);
               u2='y'+str(zz);
               a1=0;a2=0;
               uu=np.array([u1,u2]);
               rec=np.array([a1,a2]);
               usize=usize+2;               #Use to Build Assembly matrix
               u.append(uu);
               aDOF.append(rec);
               globnodes1=0;
               #globnodesextra.append(globnodes1);
        
       
        globnodesextra.append(globnodes1);
        globnodes.append(globnode);
        if rec[0] == 0:  
            globU.append(uu);
        else:
            tempU.append(uu);
            tempU.append(rec);
            tempU2=np.array(tempU);
            tempU1=tempU2.flatten();
            globU.append(tempU1);
            #globU.append(uu)
            #globU.append(rec);
            
            
        
        zz=zz+1;
        
  return globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF