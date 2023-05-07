import numpy as np


def Totbodyforcefunc(i,Bmat,Bmatenr,element,elementnoc,wi,t,det,bodforce,totbodyforce,N1,N2,
                     N3,N4,c1,c2,xi1,xi2,phivalues):
     sizoff=np.shape(Bmat);
     sizoff1=np.shape(Bmatenr);
    # for i in range(len(element)): 
     #fbstd=np.zeros(sizoff[1]);
         #fbenr=np.zeros(sizoff1[1]);  
     def bodyforce(i,elemtyp,j):
        N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j]);
        if elemtyp=="Std":
           Nstd=np.array([[N11,0,N12,0,N13,0,N14,0],[0,N11,0,N12,0,N13,0,N14]]);
           return Nstd

        elif elemtyp=="Enr" :
           phi=phivalues[i];
           Nenr=np.array([[N11,0,N12,0,N13,0,N14,0,(N11*phi[0]),0,(N12*phi[1]),0,\
                           (N13*phi[2]),0,(N14*phi[3]),0],[0,N11,0,N12,0,N13,0,N14,0,(N11*phi[0]),0,\
                                           (N12*phi[1]),0,(N13*phi[2]),0,(N14*phi[3])]]);
           return Nenr         
     if i in elementnoc:
            fbstd=np.zeros(sizoff[1]);
            elemtyp="Std"
            for j in range(4):
                  fb1=np.float32(wi[j]*np.matmul(np.transpose(bodyforce(i,elemtyp,j)),bodforce)*t*det);
                  fbstd=fb1+fbstd;
                  fb=fbstd.reshape(-1,1,order='F');
        
        
    
     else :
            fbenr=np.zeros(sizoff1[1]);
            for j in range(4):
                elemtyp="Enr"
                fb1=np.float32(wi[j]*np.matmul(np.transpose(bodyforce(i,elemtyp,j)),bodforce)*t*det);
                
                fbenr=fb1+fbenr;
                fb=fbenr.reshape(-1,1,order='F');
            
            
         
     return fb