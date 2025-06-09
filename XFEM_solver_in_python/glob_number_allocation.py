import numpy as np

globnumber=[];elementswglobnode=[];
def glob_number(noc,nor,count,elementc,u,aDOF) : 
 for i in range(noc):  
    for j in range(nor):
        #RR=nor+1;
        RR=nor+1;
         #CC=noc+1;
       # CC=noc+1;
        check=[];
        for k in range(len(elementc)):  #To check the new(variable) list,we are using for loop(k):
            check1=np.array_equal(count,elementc[k]);
            check.append(check1);
            if check1==True:
                break;
        if any(check) == False:
              gnode1=u[j+((i+1)*RR)];
              a=j+((i+1)*RR);
              gnode2=u[j+(i*(RR))];
              b=j+(i*(RR));
              gnode3=u[(j+1)+(i*(RR))];
              c=(j+1)+(i*(RR));
              gnode4=u[(j+1)+(RR*(i+1))];
              d=(j+1)+(RR*(i+1));
              globnumber1=np.array([a,b,c,d]);
              elementswglobnode1=np.array([gnode1,gnode2,gnode3,gnode4]);
              elementswglobnode.append(elementswglobnode1);
              globnumber.append(globnumber1);
             
              
        else :
             gnode1=u[j+((i+1)*RR)];
             a=j+((i+1)*RR);
             gnode2=u[j+(i*(RR))];
             b=j+(i*(RR));
             gnode3=u[(j+1)+(i*(RR))];
             c=(j+1)+(i*(RR));
             gnode4=u[(j+1)+(RR*(i+1))];
             d=(j+1)+(RR*(i+1));
             enrnode1=aDOF[j+(RR)+(i*(RR))];
             enrnode2=aDOF[j+(i*(RR))];
             enrnode3=aDOF[(j+1)+(i*(RR))];
             enrnode4=aDOF[(j+1)+(RR)+(i*(RR))];
             globnumber2=np.array([a,b,c,d]);
             elementswglobnode2=np.array([gnode1,enrnode1,gnode2,enrnode2,gnode3,enrnode3,gnode4,enrnode4]);
             elementswglobnode.append(elementswglobnode2);
             globnumber.append(globnumber2);
             
        
        
        count=count+1;
 return elementswglobnode,globnumber
