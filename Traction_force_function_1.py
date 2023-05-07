import numpy as np  ##Dupicate copy


def Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,wi,belement,
              N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues,Edge,cornerelements):
   sizoff=np.shape(Bmatstd);
   sizoff1=np.shape(Bmatenr);
   #for i in range(len(element)): 
   fbstd=np.zeros(sizoff[1]);        #This is for elements without forces(standard)
   fbenr=np.zeros([3,16]);
   #fbenr=np.zeros(sizoff1[1]);        #This is for elements without forces(Enriched)
       #fbenr=np.zeros(sizoff1[1]); 
   def shapefunctract(i,elemtyp,j):
        N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j]);
        N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j]);
        if elemtyp=="Std":
           Nstd=np.array([[N11,0,N12,0,N13,0,N14,0],[0,N11,0,N12,0,N13,0,N14]]);
           return Nstd

        elif elemtyp=="Enr" :
           phi=phivalues[i];
           Nenr=np.array([[N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2]),0,N14,0,(N14*phi[3]),0],\
                           [0,N11,0,(N11*phi[0]),0,N12,0,(N12*phi[1]),0,N13,0,(N13*phi[2]),0,N14,0,\
                                          (N14*phi[3])]]);
           return Nenr             
   if i in elementnoc and i in belement:
           elemtyp="Std"
           #fbstd1=np.zeros(sizoff[1]);
           for j in range(4):
              N=shapefunctract(i, elemtyp, j)
              dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])  #Use this method(correct)
              dNdxi1=np.round(dNdxi1,2)
              kk=element[i]
              Jac1=np.float32(np.matmul(dNdxi1,kk));           #jacobi matrix
              det=np.linalg.det(Jac1);   
              fb1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*det);
              fbstd=fb1+fbstd;
           if Edge == 'Edge1':
               # if i !=cornerelements[0] and i !=cornerelements[1]:
                   fbstd[0:2]=0;
                   fbstd[6:8]=0;
                   #fbstd[4:8]= fbstd[4:8]
          
           elif Edge == 'Edge2':
                 #if i!=cornerelements[1] and i!=cornerelements[3] :
                     fbstd[0:4]=0;
                
           elif Edge == 'Edge3':
                   # if i!=cornerelements[2] and i!=cornerelements[3] :
                        fbstd[2:6]=0;
                    #    fbstd[0:2]=fbstd[0:2]/2
                    #    fbstd[6:8]=fbstd[6:8]/2
                  
           elif Edge == 'Edge4':
               #     if i!=cornerelements[0] and i!=cornerelements[2] :
                        if i==0:
                            fbstd[0:8]=0;
                        else:    
                            fbstd[4:8]=0;
               #         fbstd[0:4]=fbstd[0:4]/2;
               
                    
           #fbstd1=fb1+fbstd1;
           fb=fbstd.reshape(-1,1,order='F');
        
        
    
   elif i in elementc and i in belement :     ###Pending
           #fbenr1=np.zeros(sizoff1[1]);
           elemtyp="Enr"
           for j in range(4):
               #elemtyp="Enr"
               N=shapefunctract(i, elemtyp, j)
               dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])  #Use this method(correct)
               dNdxi1=np.round(dNdxi1,2)
               kk=element[i]
               Jac1=np.float32(np.matmul(dNdxi1,kk));           #jacobi matrix
               det=np.linalg.det(Jac1);   
               fb1=np.float32(wi[j]*np.matmul(np.transpose(N),tractforces)*t*det);
               fbenr=fb1+fbenr;
           if Edge == 'Edge1':
                     fbenr[4:12]=0;
           elif Edge == 'Edge2':
                     fbenr[0:8]=0;
           elif Edge == 'Edge3':
                     fbenr[0:4]=0;
                     fbenr[12:16]=0;
           else:
                     fbenr[8:16]=0;
                     
                     
           #fbenr1=fb1+fbenr1;
           fb=fbenr.reshape(-1,1,order='F');
               
   else:
        if i in elementnoc and i not in belement:
            fb=fbstd.reshape(-1,1,order='F');
        elif i in elementc and i not in belement:
            fb=fbenr.reshape(-1,1,order='F');
            
            
            
   #totTractforce.append(fb);
   return fb
