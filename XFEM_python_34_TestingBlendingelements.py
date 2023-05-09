####  Author : Nandha Gopal Mariappan######
####  Title  : Programming an XFEM module for investigating the holes and inclusions in a 2D plate###
####  Date : 16-11-2022     ###############

import math
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from Body_force_function import Totbodyforcefunc
from Body_force_assembly import Totbodyforceassm
from Boundary_elements import Tractionforce
from Traction_force_function_1 import Tractforce  ###Testing
from Traction_force_assembly import Tottractforceassm
from Boundary_nodes import fixednodes
from Displacement_stress_Strain import Disp_stress_strain
from Dirichlet_Boundarycondition import DirichletBC
#from Determinant_function import my_funcdet
from LevelSetFunction import LevelSet
from Rectangular_sub_grid_nodes_1 import Partitioning_method
#from Determinant_function_partitioning import my_funcdetforpartition
from Building_elements_1 import creatingelement
from B_matrix_for_std import my_func
from B_matrix_for_enhanced import my_funcenh
from Assembly_Stiffness_mat import assembly_K
from Naming_the_umatrix import naming_umatrix
from Stress_strain_components import stress_strain_component

L=2;                    #Length of the plate(in cm)(---->Vertical length)-->Corresponds to columns-->y=direction
br=2;                    #Breadth of the plate(in cm)(Horizontal length)-->Corresponds to rows-->x-direction
l = 0.05;                   #element size-square (in cm)
nor= int((br/l));         #Number of elements in a meshed plate in X-direction
noc= int((L/l));         #Number of elements in a meshed plate in Y-direction
norforpartition = 3
nocforpartition = 3
element=[];              #Coordinates of the nodes(origin)
elementnoc=[];             #Containing elements that is not cut by 
elementc=[];
phivalues=[];            #All phivalues of an element(4nodes) will be saved in this list
phivalforpart=[];
a1=1;b1=1;               #Coordinates of the center of the circle
r1=0.4;                    #Radius of the circle
a2=0;b2=0;                  ##a respect to breadth(b)....b respects to Length(L)
r2=0;
#DOF=[];
DOF=2;
u=[];
zz=[];
N=[];
pp=[];
globnodes=[];
elementswglobnode=[];
aDOF=[];
Kloc=[];
globU=[];          
E1=1*(10**5)         #Young's modulus in Kg/cm2(Steel plate)
E2=0.01*E1;
#E2=1*(10**5)           #Young's modulus in Kg/cm2(Hole) #To avoid singular stiffness matrices1(E2/E1=0.01)
#E2=0; 
E3=4.19*(10**5)    
nu=0.3;
globnodesextra=[];
globnumber=[];
t=1;                     #Thickness of the plate(cm)
bodforce=np.array([0,0]);  #In newton   (changed the direction and symbol)
#tractforces=np.array([39.22,0]); #In newton  (changed the direction and symbol)
totbodyforce=[];
#belement1=[];
totTractforce=[];
#F3=[];
cornerelements=[];
fixedDOFSlist1=[];
l3=[];
sigma=[];
umatrix=[];
displacement=[];
Bmatrix=[];
elementinsidehole=[];
elementinsideinc=[];
elementcuthole=[];
elementcutinc=[];
stress=[];
strain=[];
#x=[];
parlist=[];
Edges=[];
#phivaluesforinc=[];
#elementsbyinc=[];

#cornernodes=np.array([(0,0),(nor-1,0),(0,noc-1),(nor-1,noc-1)])
###To build the elements and numbering the nodes###(local nodes)
iteration=0;
for i in range(noc):  
    for j in range(nor):
        elementcall,cornerelements1=creatingelement(i,j,l,noc,nor)
        element.append(elementcall)
        if cornerelements1==1:
            cornerelements.append(iteration)
        iteration=iteration+1;
        
        
 
    
#Totelements=len(element);
    
###To check the level-set function for all nodes
for i in range(len(element)):
    temp=element[i];
    if r1>0 and r2==0:
      aa=LevelSet(temp,a1,b1,r1)
      tt=[0,0,0,0]
    else:
      tt=LevelSet(temp,a2,b2,r2)
        
    #if np.sum(a)== 4 or np.sum(a)== 0  :
    #if ([0,-1] in aa,tt) or (-1 in aa,tt):
    if np.sum(aa)==4 or np.sum(tt)==4:    ##Changed 'and' to 'or' because of removal of inclusion
        elementnoc.append(i);
    else:
      if np.sum(aa)<4 or np.sum(tt)<4:
        if all(phi==-1 for phi in aa)== True:
            elementnoc.append(i);
            elementinsidehole.append(i)
        elif all(phi==-1 for phi in tt)== True:
            elementnoc.append(i);
            elementinsideinc.append(i)
        elif 0 in aa and -1 not in aa :#or 0 in tt and -1 not in tt:
            elementnoc.append(i)
        else:
            elementc.append(i);
            if -1 in aa:
                elementcuthole.append(i);
            elif -1 in tt:
                elementcutinc.append(i);
    #else:
        #elementnoc.append(i);
    if -1 in aa or 0 in aa:
       phivalues.append(aa);                    #This phivalues variable will tell about the nodes which lies 
   # elif -1 in tt or 0 in tt:                              #in an element whether it is inside or outside the curve
    #   phivalues.append(tt);
    else:
       phivalues.append(aa)
                      
                                          
########For partioning the element with curves for better accuracy#############
for i in range(len(elementc)):
     parlist1=Partitioning_method(i,l,elementc,element) 
     parlist.append(parlist1)
                             
####For partitioned elements,calculating the phi values#######################
for i in range(len(elementc)):
     phivalforpart1=[];
     bb=elementc[i];
     subgrids=parlist[i]
     if bb in elementcuthole:
        for j in range(len(subgrids)):
           b3=LevelSet(subgrids[j],a1,b1,r1)
           phivalforpart1.append(b3)
        phivalforpart.append(phivalforpart1)
     elif bb in elementcutinc:
        for j in range(len(subgrids)):
           b3=LevelSet(subgrids[j],a2,b2,r2)
           phivalforpart1.append(b3)
        phivalforpart.append(phivalforpart1)
                                              
                                            
                                            
                                            
  
for i in range(len(elementc)):     #To remove the duplicate coordinates from the enriched elements
     for j in range(4):            #which helps to differentiate the enriched global nodes from the                          
       nn=(element[elementc[i]]);  #ordinary global nodes
       cp=np.array(nn[j]);
       pp.append(cp);
       
new = np.unique(pp,axis=0);
#print (new);

zz=0;
aDOFsiz=0;
usize=0;

globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF=naming_umatrix(noc,nor,l,new,zz,aDOFsiz,usize,globnodesextra,globnodes,globU,u,aDOF)
#globnode,aDOFsiz,usize,globnodesextra,globU=naming_umatrix(noc,nor,l,new,zz,aDOFsiz,usize):
      
#globnodes=np.array(globnodes2)
flatten1 = list(np.concatenate(globU).flat);       #Flatten1 has all the 2d array to a single list

#Using the above function,it helps to break the 2d array to an 1d array
flatten = np.array([flatten1]).T;                  #Transpose of an flatten1 matrix(contains
#information of globU)


count=0;        
for i in range(noc):  
    for j in range(nor):
        #RR=nor+1;
        RR=nor+1;
        #CC=noc+1;
        CC=noc+1;
        check=[];
        for k in range(len(elementc)):  #To check the new(variable) list,we are using for loop(k):
            check1=np.array_equal(count,elementc[k]);
            check.append(check1);
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
        
        

##############2D Gaussian quadrature with points and weights#################
xi1 = np.array([(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3),math.sqrt(1/3)]);
xi2 = np.array([(-math.sqrt(1/3)),(math.sqrt(1/3)),(-math.sqrt(1/3)),(math.sqrt(1/3))]);
wi = np.array([1,1,1,1]);

#Shape Functions:-
c1 = sp.Symbol('c1');
c2 = sp.Symbol('c2');

N1=((1-c1)*(1-c2))/4;
N2=((1+c1)*(1-c2))/4;      
N3=((1+c1)*(1+c2))/4;
N4=((1-c1)*(1+c2))/4;

Nstd=np.array([[N1,0,N2,0,N3,0,N4,0],[0,N1,0,N2,0,N3,0,N4]])

'''
def my_func(i,j):            #To build the B-matrix for standard elements
        
         cc=i;dd=j
         dNdxi,detj,Jac = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
         #dNdxi,detj,Jac = my_funcdet(cc,dd);
         #detj=np.linalg.det(Jac)                 #To find determinant of the jacobi matrix
         invj=np.linalg.inv(Jac)                 #To find inverse of the jacobi matrix
         #c2=xi2[i];
         Bstd1=np.matmul(invj,dNdxi);
         
         Bstd =np.array([[Bstd1[0,0],0,Bstd1[0,1],0,Bstd1[0,2],0,Bstd1[0,3],0],\
                        [0,Bstd1[1,0],0,Bstd1[1,1],0,Bstd1[1,2],0,Bstd1[1,3]],\
                           [Bstd1[1,0],Bstd1[0,0],Bstd1[1,1],Bstd1[0,1],\
                             Bstd1[1,2],Bstd1[0,2],Bstd1[1,3],Bstd1[0,3]]]);
             
         #Bstd=np.matmul(invj,Bstd1);
         
         
         return Bstd,detj
''' 
'''    
def my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,Dvalues1):          #To build the B-matrix for enriched elements
        
         cc=i;dd=j
         #dNdxi,detj,Jac = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element);  (Last Update:27/03/2023)-line
         detjenr = my_funcdet(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
         dNdxi,detjenrpart,Jac = my_funcdetforpartition(cc,dd,xi1,xi2,c1,c2,N1,N2,N3,N4,element,parlist,part,k);
         
       
         invj=np.linalg.inv(Jac)                #To find inverse of the jacobi matrix
         #print (detj)
         
         phivalues2=phivalforpart[part]
         phivalues1=phivalues2[k]
         ##New addition(01-04-23)
         if phivalues1[j]==1:
             D=Dstd
             DD="std"
         elif phivalues1[j]==0:
             if i in elementcuthole:
                    D=D1;
                    DD="Dhole"
             else:
                    D=D2;
                    DD="Dinc"
         else:
             if i in elementcuthole:
                    D=D1;
                    DD="Dhole"
             else:
                    D=D2;
                    DD="Dinc"
             
         Benr1=np.matmul(invj,dNdxi);
         Dvalues1.append(DD);
         
        
          
          
         Benr=np.array([[(Benr1[0,0]),0,Benr1[0,0]*float(phivalues1[0]),0,Benr1[0,1],0,Benr1[0,1]*float(phivalues1[1]),0,Benr1[0,2],0,\
                      Benr1[0,2]*float(phivalues1[2]),0,Benr1[0,3],0,Benr1[0,3]*float(phivalues1[3]),0],\
                      [0,Benr1[1,0],0,Benr1[1,0]*float(phivalues1[0]),0,Benr1[1,1],0,Benr1[1,1]*float(phivalues1[1]),0,Benr1[1,2],0,Benr1[1,2]*float(phivalues1[2]),\
                       0,Benr1[1,3],0,Benr1[1,3]*float(phivalues1[3])],\
                         [Benr1[1,0],Benr1[0,0],Benr1[1,0]*float(phivalues1[0]),Benr1[0,0]*float(phivalues1[0]),Benr1[1,1],Benr1[0,1],\
                          Benr1[1,1]*float(phivalues1[1]),Benr1[0,1]*float(phivalues1[1]),Benr1[1,2],Benr1[0,2],Benr1[1,2]*float(phivalues1[2]),Benr1[0,2]*float(phivalues1[2]),\
                          Benr1[1,3],Benr1[0,3],Benr1[1,3]*float(phivalues1[3]),Benr1[0,3]*float(phivalues1[3])]]);
         return Benr,detjenr,detjenrpart,D,Dvalues1;
'''     
#Kenr=0;
#Kstd=0;
part=0;  #Counter for partition elements
Dvalues=[];
Dvalues2=[];

C=E1/((1-nu)*(1-(2*nu)))                  #Plane strain condition
#C=(E1/1-(nu**2));
#Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate
Dstd=np.array([[C*(1-nu),(C*nu),0],[(C*nu),C*(1-nu),0],[0,0,(((1-(2*nu))/2)*C)]]);     #Plane strain condition

C1=E1/((1-nu)*(1-(2*nu)))                  #Plane strain condition
#C1=(E2/1-(nu**2));
#D1=np.array([[C1,(C1*nu),0],[(C1*nu),C1,0],[0,0,(((1-nu)/2)*C1)]]);  #Material Tangent Matrix for a hole
D1=np.array([[C1*(1-nu),(C1*nu),0],[(C1*nu),C1*(1-nu),0],[0,0,(((1-(2*nu))/2)*C1)]]);     #Plane strain condition

C2=(E3/1-(nu**2));
D2=np.array([[C2,(C2*nu),0],[(C2*nu),C2,0],[0,0,(((1-nu)/2)*C2)]]);   #Material Tangent matrix for an inclusion

#Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);
#Denr=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);
totDOF=usize+aDOFsiz;
#K=np.zeros((totDOF,totDOF));      
for i in range(len(element)):        #Pending-- weights need to be implemented in stiffness matrix - Updated
      Kenr=0;
      Kstd=0;
      sigmaenr=0;
      sigmastd=0;
      sigmaenrfl=0;
      if i in elementc:    ##To create enriched B-matrix(Updated one)
        Dvalues3=[];
        for k in range(len(parlist[part])):
            Dvalues1=[];
            for j in range(4):
                  #Bmatenr,detstd,D= my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc,phivalforpart,phivalues); ##Old(01-04-23)
                  Bmatenr,detstd,detenrpart,D,Dvalues2= my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,Dvalues1,
                                                                    xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc,phivalforpart);  ##original
                  KK2=np.matmul(np.transpose(Bmatenr),D);
                  sigmaenr1=np.matmul(D,Bmatenr)*detstd
                  K2=wi[j]*np.matmul(KK2,Bmatenr)*detenrpart*detstd;    ###Detstd gives determinant for the single element
                  #K2=wi[j]*np.matmul(KK2,Bmatenr)*detstd; 
                  Kenr=K2+Kenr;
                  sigmaenr=sigmaenr1+sigmaenr;
                  K4=(Kenr);
            Dvalues3.append(Dvalues2)
            sigmaenrfl=(sigmaenr)+sigmaenrfl;                            ###multiplication of D and Bmatenr for all single subgrids
        Bmatrix.append(Bmatenr);
        sigma.append(sigmaenrfl);
        Dvalues.append(Dvalues3);
        part=part+1;
           
      else:
            for j in range(4): ##To create standard B-matrix
               #J,Bmat= my_func(i,j);
               #Jac=np.matmul(J,element[i]);
               Bmatstd,detstd= my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
               if i in elementinsidehole:##New line
                   KK2=np.matmul(np.transpose(Bmatstd),D1);
                   sigmastd1=np.matmul(D1,Bmatstd)
               elif i in elementinsideinc:
                   KK2=np.matmul(np.transpose(Bmatstd),D2);
                   sigmastd1=np.matmul(D2,Bmatstd)
               else:   
                   KK2=np.matmul(np.transpose(Bmatstd),Dstd);
                   sigmastd1=np.matmul(Dstd,Bmatstd)
               K1=wi[j]*np.matmul(KK2,Bmatstd)*detstd;
               Kstd=K1+Kstd;
               sigmastd=sigmastd1+sigmastd;
               K4=(Kstd);
            Bmatrix.append(Bmatstd); 
            sigma.append(sigmastd); 
              
           
      Kloc.append((K4));
      #sigma.append(KK2);  ###Sigma calc is wrong(Detected :12-04-23)
      det=detstd;
      



##################Assembly of K-matrix(Stiffness matrix)#######################
K=np.zeros((totDOF,totDOF));
K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)
'''
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
'''           
flatteig,flatt1=np.linalg.eig(Kloc[23]);#To check the whether the stiffness matrix is singular or not(it is a singular matrix)
flatteig1,flatt2=np.linalg.eig(Kloc[26]);
sum1=np.sum(flatteig)
sum2=np.sum(flatteig1)
#Flatt=np.linalg.solve(K,F)

##########################To build the element body force matrix:-###################################
##Body forces acting on the plate will be same for all the elements because forces acts equally
##on all elements under gravity.

for i in range(len(element)):
   # Shapefunc=Shapefunction(N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues)
    Totalforcebod = Totbodyforcefunc(i,Bmatstd,Bmatenr,element,elementnoc,wi,t,det,bodforce,totbodyforce,
                                     N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues)
    totbodyforce.append(Totalforcebod);
##Assembly of F matrix(body force)
Fbod1=np.zeros(np.shape(flatten));
for i in range(len(element)):
    fbloc1=totbodyforce[i];
    Fbod=Totbodyforceassm(i,globnumber,elementnoc,DOF,globnodesextra,fbloc1,Fbod1)
          
#########################   Traction force         ###########################################
           
###Ask from the user for applying Traction force regarding which edge needed for them
Ftract1=np.zeros(np.shape(flatten));

Answer=input("Need to Apply Traction force,please type Yes or no :")
if Answer == 'Yes':
 while True:
  #Answer=input("Need to Apply Traction force,please type Yes or no :")
  belement1=[];
  if Answer == 'No':
      belement=[];
      #F3=0;
      break;
  else:
     Edge=input("Enter the edge number for applying Traction force : ")  
     Tractforcex=input("Enter the traction force in the x-direction : ")
     Tractforcey=input("Enter the traction force in the y-direction : ")
     tractforces=np.array([float(Tractforcex),float(Tractforcey)]);
#Edge="Edge1"   #Temporary
  
     if Edge == "Edge1": 
          belement1.append(Tractionforce(nor,noc,"Edge1",globnodes,globnumber,l))
     elif Edge == "Edge2":
          belement1.append(Tractionforce(noc,nor,"Edge2",globnodes,globnumber,l))
     elif Edge == "Edge3":
          belement1.append(Tractionforce(nor,noc,"Edge3",globnodes,globnumber,l))
     elif Edge == "Edge4":
          belement1.append(Tractionforce(noc,nor,"Edge4",globnodes,globnumber,l))
     else:
        print("The entered edge has not been found amoung those edges")
     
  #belement2.append(belement1)
     F3=[];
     belement=list(np.concatenate(belement1).flat)
     for i in range(len(element)):
            F2=Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,wi,belement,
                      N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues,Edge,cornerelements)
            F3.append(F2);
            
     for i in range(len(element)):
      fbloc2=F3[i];
      Ftract=Tottractforceassm(i,globnumber,elementnoc,DOF,globnodesextra,fbloc2,Ftract1)  
   
  ii = input("Want to give Traction force for another edge Type(Yes or No): ")   #Yes or no
  if ii == "No":
    break; 
'''
elif Answer == 'No':
    belement=[];
    F3=0;
    '''
####To find corner elements for traction forces####
#cornerelements=np.array([0,9,90,99])

###Put the Tractforce inside the belement#####494 to 500
'''
if Answer == 'Yes':
    for i in range(len(element)):
           F2=Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,det,wi,belement,
                     N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues,Edge,cornerelements)
           F3.append(F2);
else:
       F3=0       
###Put the Tractforce inside the belement#####494 to 500 
'''   
  #Ftract1=np.zeros(np.shape(flatten));
'''
if F3==0:
    Ftract=np.zeros(np.shape(flatten));
else:
    for i in range(len(element)):
      fbloc2=F3[i];
      Ftract=Tottractforceassm(i,globnumber,elementnoc,DOF,globnodesextra,fbloc2,Ftract1)           
'''   
#  ii = input("Want to give Traction force for another edge : ")   #Yes or no
 # if ii == "no":
  #  break;         
         
F=Fbod+Ftract;


######For finding the fixed nodes(global number)########which helps in Dirichllet BC
#fixed_edge=input("Enter the name of the fixed edge : ")
BCnodes=[]; ##Checking for boundary nodes
Constraints=input("Are you need to fix/constraints the Edges? Type:(Yes or No)")
if Constraints == 'No':
    fixedDOFSlist1=[]
else:
 while True: 
  fixed_edge=input("Enter the name of the fixed edge : ")
  boundarynodes=[]
  #fixed_edge="Edge3"     #Temporary
  #if fixed_edge =="No":    #Temp
     # break;
  if fixed_edge == "Edge1": 
         bnodes = fixednodes(nor,noc,"Edge1",globnodes,globnumber,l)
  elif fixed_edge == "Edge2":
         bnodes = fixednodes(noc,nor,"Edge2",globnodes,globnumber,l)
  elif fixed_edge == "Edge3":
         bnodes = fixednodes(nor,noc,"Edge3",globnodes,globnumber,l)
  elif fixed_edge == "Edge4":
         bnodes = fixednodes(noc,nor,"Edge4",globnodes,globnumber,l)
  else:
      print("The entered edge has not been found amoung those edges")       
      
  boundarynodes.extend(bnodes)   ##Bnodes and boundary nodes are same
  BCnodes.append(boundarynodes)         ##Checkpoint : BCnodes
  
  #################Applying Dirchelet Boundary Conditions###################
  FixedDOFS3=input("Do you need to constrain the nodes of respective edge in the x-dir or y-dir or both?")
  for i in range(len(boundarynodes)):
     j=boundarynodes[i]
  
     if FixedDOFS3 == "x":
         fixedDOFS2=globU[j]
         if fixedDOFS2.size == 2:
             fixedDOFS4=np.delete(fixedDOFS2,1)   #since x is fixed,remove the y-values
         else:
             fixedDOFS4=np.delete(fixedDOFS2,[1,3])
             #fixedDOFS4=list(np.concatenate(fixedDOFS4).flat)
     elif FixedDOFS3 == "y":
         fixedDOFS2=globU[j]
         if fixedDOFS2.size == 2:
             fixedDOFS4=np.delete(fixedDOFS2,0)   #since y is fixed,remove the x-values
         else:
             fixedDOFS4=np.delete(fixedDOFS2,[0,2])
             #fixedDOFS4=list(np.concatenate(fixedDOFS4).flat)
     elif FixedDOFS3 == "both":
         fixedDOFS4=globU[j]
     
     fixedDOFSlist1.append(fixedDOFS4)
     fixedDOFSlist=list(np.concatenate(fixedDOFSlist1).flat)
     
  for i in range(len(fixedDOFSlist)):
        fixDOF=fixedDOFSlist[i];
        K,l3=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3)
        
  ii = input("Want to fix another edge,Please Type Yes or No : ")   #Yes or no
  if ii == "No":
      break;

'''
if fixed_edge == 'No':
    fixedDOFSlist=[];
else:           
    for i in range(len(fixedDOFSlist)):
       fixDOF=fixedDOFSlist[i];
       K,l3=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3)
'''

fixedDOFSlist=['x'+str(0),'y'+str(0),'x'+str(1560),'y'+str(1560)]
for i in range(len(fixedDOFSlist)):
      fixDOF=fixedDOFSlist[i];
      K,l3=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3)

#BCnodes1=list(np.concatenate(BCnodes).flat)    ###Checkpoint for boundarynodes..BCnodes 1should match with fixedDOFlist
  #ii = input("Want to fix another edge : ")   #Yes or no
  #if ii == "no":
  #    break;
      

  
flatten=np.linalg.solve(K,F);  


###To calculate displacements,stresses,strains for all elements:-
increment=0;
for i in range(len(element)):
    disp=umatrix[i];   
    disp3,actdisp=Disp_stress_strain(i,disp,flatten,sigma,Bmatrix,elementc,elementnoc,elementinsidehole,phivalues)         ##Stress and strain will come from this block
    displacement.append(actdisp);
    stress1,strain1=stress_strain_component(i,disp,flatten,element,elementc,elementnoc,elementinsidehole,Dstd,D1,phivalues,\
                                            phivalforpart,increment)
    stress.append(stress1)
    strain.append(strain1)
 

 
xvalues2=np.linspace(0,br,nor)   ###Arranging elements
yvalues2=np.linspace(0,L,noc)    ###Arranging elements
xv1, yv1 = np.meshgrid(xvalues2, yvalues2)
zvalues2=np.zeros([nor,noc])
zvalues3=np.zeros([nor,noc])
loop=0;
for i in range(len(xvalues2)):
    for j in range(len(yvalues2)):
        stressforeach=stress[loop]
        zvalues2[i][j]=stressforeach[0]
        zvalues3[i][j]=stressforeach[1]
        loop=loop+1;
      

#plt.plot(xv, yv, marker='.', color='k', linestyle='none')
fig,ax2=plt.subplots(1,1)
fig,ax3=plt.subplots(1,1)
ax2.set_title('Stress contour in x-direction (σx) (kg/cm2),',fontweight='bold')
ax3.set_title('Stress contour in y-direction (σy) (kg/cm2),',fontweight='bold')
ax2.set_xlabel('Breadth of the plate',fontweight='bold')
ax2.set_ylabel('Length of the plate',fontweight='bold')
ax3.set_xlabel('Breadth of the plate',fontweight='bold')
ax3.set_ylabel('Length of the plate',fontweight='bold')
cp2 = ax2.contourf(xv1, yv1, zvalues2,cmap='jet') 
cp3 = ax3.contourf(xv1, yv1, zvalues3,cmap='jet')
fig.colorbar(cp2)
fig.colorbar(cp3)
plt.autoscale(False)
plt.show()
plt.show() 


dispcontx=0
dispconty= 0  
xvalues=np.linspace(0,br,nor)   ###Arranging elements
yvalues=np.linspace(0,L,noc)    ###Arranging elements
xv, yv = np.meshgrid(xvalues, yvalues)
zvalues=np.zeros([nor,noc])
zvalues1=np.zeros([nor,noc])
loop1=0;
for i in range(len(xvalues)):
    for j in range(len(yvalues)):
        dispall=displacement[loop1]
        zvalues[i][j]=dispall[0]
        zvalues1[i][j]=dispall[1]
        loop1=loop1+1;
        

#plt.plot(xv, yv, marker='.', color='k', linestyle='none')
fig,ax=plt.subplots(1,1)
fig,ax1=plt.subplots(1,1)
#ax.set_title('Contour plot of x-directional displacement (in mm)',fontweight='bold')
#ax1.set_title('Contour plot of y-directional displacement (in mm)',fontweight='bold')
ax.set_xlabel('Breadth of the plate',fontweight='bold')
ax.set_ylabel('Length of the plate',fontweight='bold')
ax1.set_xlabel('Breadth of the plate',fontweight='bold')
ax1.set_ylabel('Length of the plate',fontweight='bold')
cp = ax.contourf(xv, yv, zvalues,cmap='jet') ##Problem in contour plot
cp1 =ax1.contourf(xv, yv, zvalues1,cmap='jet') ##Problem in contour plot
fig.colorbar(cp)
fig.colorbar(cp1)
plt.autoscale(False)
plt.show()
plt.show()       
