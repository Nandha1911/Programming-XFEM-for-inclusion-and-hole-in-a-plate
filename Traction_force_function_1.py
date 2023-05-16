import numpy as np  ##Dupicate copy


def Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,wi,belement,
              N1,N2,N3,N4,c1,c2,phivalues,Edge,l,F2):
 
 for i in range(len(element)):
   sizoff=np.shape(Bmatstd);
   sizoff1=np.shape(Bmatenr);
   #for i in range(len(element)): 
   c=int((sizoff[1])/2)
   d=int((sizoff1[1])/2)
   ftstd1=np.zeros(c);        #This is for elements without forces(standard)
   ftenr1=np.zeros(d);
   ftstd=np.zeros((sizoff[1]));        #This is for elements without forces(standard)
   ftenr=np.zeros((sizoff1[1]));
  # xi1 = np.array([-1,1,1,-1]);
  # xi2 = np.array([-1,-1,1,1]);
   #fbenr=np.zeros(sizoff1[1]);        #This is for elements without forces(Enriched)
       #fbenr=np.zeros(sizoff1[1]); 
   def shapefunctract(i,elemtyp,j,Edge1):
        if Edge1=='Edge1':
            xi1 = np.array([1,1]);
            xi2 = np.array([-1,1]);
            N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j]);
            N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j]);
            if elemtyp=="Std":
               Nstd=np.array([[N12,0,N13,0],[0,N12,0,N13]]);
               return Nstd

            elif elemtyp=="Enr" :
               phi=phivalues[i];
               Nenr=np.array([[N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2]),0],\
                              [0,N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2])]]);
               return Nenr
        elif Edge1 == 'Edge2':
            xi1 = np.array([1,-1]);
            xi2 = np.array([1,1]);
            N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j]);
            N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j]);
            if elemtyp=="Std":
               Nstd=np.array([[N13,0,N14,0],[0,N13,0,N14]]);
               return Nstd

            elif elemtyp=="Enr" :
               phi=phivalues[i];
               Nenr=np.array([[N13,0,(N13*phi[2]),0,N14,0,(N14*phi[3]),0],\
                              [0,N13,0,(N13*phi[2]),0,N14,0,(N14*phi[3])]]);
               return Nenr
        elif Edge1 == 'Edge3':
               xi1 = np.array([-1,-1]);
               xi2 = np.array([-1,1]);
               N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j]);
               N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j]);
               if elemtyp=="Std":
                  Nstd=np.array([[N11,0,N14,0],[0,N11,0,N14]]);
                  return Nstd

               elif elemtyp=="Enr" :
                  phi=phivalues[i];
                  Nenr=np.array([[N11,0,(N11*phi[0]),0,N14,0,(N14*phi[3]),0],\
                                 [0,N11,0,(N11*phi[0]),0,N14,0,(N14*phi[3])]]);
                  return Nenr
        elif Edge1 == 'Edge4':
            xi1 = np.array([-1,1]);
            xi2 = np.array([-1,-1]);
            N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j]);
            N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j]);
            if elemtyp=="Std":
               Nstd=np.array([[N11,0,N12,0],[0,N11,0,N12]]);
               return Nstd

            elif elemtyp=="Enr" :
               phi=phivalues[i];
               Nenr=np.array([[N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1]),0],\
                              [0,N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1])]]);
               return Nenr
           #Nstd=np.array([[N12,0,N13,0],[0,N12,0,N13]]);
         #   for j in range(2):
       # N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j]);
        #N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j]);
       # N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j]);
       # N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j]);
     
       
       
       
       # if elemtyp=="Std":
         #  Nstd=np.array([[N11,0,N12,0,N13,0,N14,0],[0,N11,0,N12,0,N13,0,N14]]);
        #   return Nstd

       # elif elemtyp=="Enr" :
      #     phi=phivalues[i];
        #   Nenr=np.array([[N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2]),0,N14,0,(N14*phi[3]),0],\
         #                  [0,N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2]),0,N14,0,\
        #                                  (N14*phi[3])]]);
        #   return Nenr             
   if i in elementnoc and i in belement:
           elemtyp="Std"
           #fbstd1=np.zeros(sizoff[1]);
           for j in range(2):
              N=shapefunctract(i, elemtyp, j,Edge)
              #dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])  #Use this method(correct)
              #dNdxi1=np.round(dNdxi1,2)
              #kk=element[i]
             # Jac1=np.float32(np.matmul(dNdxi1,kk));           #jacobi matrix
             #det=np.linalg.det(Jac1);   
              #fb1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*det);
              ft1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*(l/2));
              ftstd1=ft1+ftstd1;
           if Edge == 'Edge1':
                   #ftstd[0:2]=0;
                   #ftstd[6:8]=0;
                   ftstd[2:6]=ftstd1
          
           elif Edge == 'Edge2':
                     #ftstd[0:4]=0;
                     ftstd[4:8]=ftstd1
                
           elif Edge == 'Edge3':
                        #ftstd[2:6]=0;
                        ftstd[0:2]=ftstd1[0:2]
                        ftstd[6:8]=ftstd1[2:4]
           elif Edge == 'Edge4':
                        ftstd[0:4]=ftstd1
                        
               
               
                    
           #fbstd1=fb1+fbstd1;
           ft=ftstd.reshape(-1,1,order='F');
        
        
    
   elif i in elementc and i in belement :     ###Pending
           #fbenr1=np.zeros(sizoff1[1]);
           elemtyp="Enr"
           for j in range(2):
               #elemtyp="Enr"
               N=shapefunctract(i, elemtyp, j,Edge)
               #dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])  #Use this method(correct)
              # dNdxi1=np.round(dNdxi1,2)
              # kk=element[i]
              # Jac1=np.float32(np.matmul(dNdxi1,kk));           #jacobi matrix
              # det=np.linalg.det(Jac1);   
              # fb1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*det);
               ft1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*(l/2));
               ftenr1=ft1+ftenr1;
           if Edge == 'Edge1':
                     #ftenr[4:12]=0;
                     ftenr[4:12]=ftenr1;
           elif Edge == 'Edge2':
                     ftenr[8:16]=ftenr1;
           elif Edge == 'Edge3':
                     ftenr[0:4]=ftenr1[0:4];
                     ftenr[12:16]=ftenr1[4:8];
           else:
                     ftenr[0:8]=ftenr1;
                     
                     
           #fbenr1=fb1+fbenr1;
           ft=ftenr.reshape(-1,1,order='F');
               
   else:
        if i in elementnoc and i not in belement:
            ft=ftstd.reshape(-1,1,order='F');
        elif i in elementc and i not in belement:
            ft=ftenr.reshape(-1,1,order='F');
            
            
            
   F2.append(ft)
 return F2

