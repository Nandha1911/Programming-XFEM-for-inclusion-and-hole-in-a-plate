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
from Traction_force_function_2 import Tractforce  ###XFEM_36 Traction force_2 added and working correctly
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

######Inputs section###############
####___Length of the plate(in cm)(---->Vertical length)-->Corresponds to columns-->y=direction___#####
L=float(input("Enter the Length (vertical) of the plate : " ))      
#Breadth of the plate(in cm)(Horizontal length)-->Corresponds to rows-->x-direction

br=float(input("Enter the Breadth (Horizontal) of the plate : " ))                     
#Mesh size-square (in cm)

l = float(input("Enter the mesh size of the element : "))   
             
nor= int((br/l));         #Number of elements in a meshed plate in X-direction
noc= int((L/l));         #Number of elements in a meshed plate in Y-direction

#Young's modulus in Kg/cm2(Steel plate)
E=float(input("Enter the Young's modulus of the plate : "))    ##For example=1*(10**5)  

#Poisson ratio for the plate         
nu=float(input("Enter the poisson ratio for the plate : "))

element=[];              #Coordinates of the nodes(origin)
elementnoc=[];             #Containing elements that is not cut by 
elementc=[];
phivalues=[];            #All phivalues of an element(4nodes) will be saved in this list
phivalforpart=[];

#Discontinuity=input("How many discontinuities need to be added in the plate : ") ##Max:2 #Min:1

  #Coordinates of the discontinuity
a1=float(input("Enter the x coordinate of the center of the discontinuity_1 (Hole/Inclusion) : "))
b1=float(input("Enter the y coordinate of the center of the discontinuity_1(Hole/Inclusion) : "))   
            
  ##Enter the radius of the circle
r1=float(input("Enter the radius of the discontinuity_1(hole/inclusion) : "))   
 
  #Young's modulus in Kg/cm2(Steel plate)                 
E1=float(input("Enter the Young's modulus of the discontinuity_1(hole/inclusion) : ")) 

  #Poisson ratio for the discontinuity_1
nu1=float(input("Enter the poisson ratio for the discontinuity_1(hole/inclusion) : "))

Discontinuity_2=input("Do you want to provide another Discontinuity(2) in the plate (Yes or No) :")
if Discontinuity_2 == 'Yes':

#Coordinates of the center of the circle
 # a2=2.4;b2=1;                 
 # r2=0.4;
 # E2=4.19*(10**5)    
 # nu2=0.3;
  a2=float(input("Enter the x coordinate of the center of the discontinuity_2(Hole/Inclusion) : "))
  b2=float(input("Enter the y coordinate of the center of the discontinuity_2(Hole/Inclusion) : "))   
            
  ##Enter the radius of the circle
  r2=float(input("Enter the radius of the discontinuity_2(hole/inclusion) : "))   
 
  #Young's modulus in Kg/cm2(Steel plate)                 
  E2=float(input("Enter the Young's modulus of the discontinuity_2(hole/inclusion) : ")) 

  #Poisson ratio for the discontinuity_1
  nu2=float(input("Enter the poisson ratio for the discontinuity_2(hole/inclusion) : "))
  
  
condition = int((input("Do you want to apply Plane-stress or Plane-strain condition? (Type '1' for plane stress...\
                       or Type '2' for plane strain): ")))  
DOF=2;
u=[];
zz=[];
#N=[];
pp=[];
globnodes=[];
elementswglobnode=[];
aDOF=[];
Kloc=[];
globU=[];          
#E2=1*(10**5)           #Young's modulus in Kg/cm2(Hole) #To avoid singular stiffness matrices1(E2/E1=0.01)
#E2=0; 
t=1;                   #Thickness of the plate
globnodesextra=[];
globnumber=[];
bodforce=np.array([0,0]);  #In newton   (changed the direction and symbol)
#tractforces=np.array([39.22,0]); #In newton  (changed the direction and symbol)
cornerelements=[];
#fixedDOFSlist1=[];
l3=[];# Just for checking
sigma=[];
umatrix=[];
displacement=[];
#Bmatrix=[];
elementinsidehole=[];
elementinsideinc=[];
elementcuthole=[];
elementcutinc=[];
stress=[];
strain=[];
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
for i in range(len(element)):   ###Error in Level set function..Need rectification(12/05/2023)
    temp=element[i];
    aa=LevelSet(temp,a1,b1,r1)
   # tt=[];
    if Discontinuity_2 == 'Yes':
      tt=LevelSet(temp,a2,b2,r2)
    
      
    if Discontinuity_2 == 'Yes':
    ###Check for both holes
      if np.sum(aa)==4 and np.sum(tt)==4:    ##Changed 'and' to 'or' because of removal of inclusion
        elementnoc.append(i);
      else:
        if np.sum(aa)<4 or np.sum(tt)<4:
           if all(phi==-1 for phi in aa)== True:
             elementnoc.append(i);
             elementinsidehole.append(i)
           elif all(phi==-1 for phi in tt)== True:
             elementnoc.append(i);
             elementinsideinc.append(i)
           elif 0 in aa and -1 not in aa or 0 in tt and -1 not in tt:
             elementnoc.append(i)
           else:
             elementc.append(i);
             if -1 in aa:
                elementcuthole.append(i);
             elif -1 in tt:
                elementcutinc.append(i);
        #else:
          # elementnoc.append(i);
      if -1 in aa or 0 in aa:
        phivalues.append(aa);                    #This phivalues variable will tell about the nodes which lies 
      elif -1 in tt or 0 in tt:                              #in an element whether it is inside or outside the curve
        phivalues.append(tt);
      else:
        phivalues.append(aa)
    else:
      if np.sum(aa)==4:  ##Changed 'and' to 'or' because of removal of inclusion
          elementnoc.append(i);
      else:
        if np.sum(aa)<4:
          if all(phi==-1 for phi in aa)== True:
             elementnoc.append(i);
             elementinsidehole.append(i)
          elif 0 in aa and -1 not in aa :#or 0 in tt and -1 not in tt:
            elementnoc.append(i)
          else:
            elementc.append(i);
            if -1 in aa:
                elementcuthole.append(i);
      # else:
         # elementnoc.append(i);
          
      if -1 in aa or 0 in aa:
         phivalues.append(aa);                                                #in an element whether it is inside or outside the curve
      else:
         phivalues.append(aa)
    
'''
    if Discontinuity ==2:
    ###Check for both holes
      if np.sum(aa)==4 and np.sum(tt)==4:    ##Changed 'and' to 'or' because of removal of inclusion
        elementnoc.append(i);
      else:
        if np.sum(aa)<4 or np.sum(tt)<4:
           if all(phi==-1 for phi in aa)== True:
             elementnoc.append(i);
             elementinsidehole.append(i)
           elif all(phi==-1 for phi in tt)== True:
             elementnoc.append(i);
             elementinsideinc.append(i)
           elif 0 in aa and -1 not in aa or 0 in tt and -1 not in tt:
             elementnoc.append(i)
           else:
             elementc.append(i);
             if -1 in aa:
                elementcuthole.append(i);
             elif -1 in tt:
                elementcutinc.append(i);
        #else:
          # elementnoc.append(i);
      if -1 in aa or 0 in aa:
        phivalues.append(aa);                    #This phivalues variable will tell about the nodes which lies 
      elif -1 in tt or 0 in tt:                              #in an element whether it is inside or outside the curve
        phivalues.append(tt);
      else:
        phivalues.append(aa)
'''                     
                                          
########For partioning the element with curves for better accuracy#############
parlist=[];
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
    
#Kenr=0;
#Kstd=0;
part=0;  #Counter for partition elements
Dvalues=[];
Dvalues2=[];

##################___________Plane stress or plane strain conditions___________#############################
if condition == 1:  ###Plane stress condition for 1
#C=E1/((1-nu)*(1-(2*nu)))                  #Plane strain condition
  C=(E/1-(nu**2));
  Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate
#Dstd=np.array([[C*(1-nu),(C*nu),0],[(C*nu),C*(1-nu),0],[0,0,(((1-(2*nu))/2)*C)]]);     #Plane strain condition

#C1=E1/((1-nu)*(1-(2*nu)))                  #Plane strain condition
  C1=(E1/1-(nu1**2));
  D1=np.array([[C1,(C1*nu1),0],[(C1*nu1),C1,0],[0,0,(((1-nu1)/2)*C1)]]);  #Material Tangent Matrix for a hole
#D1=np.array([[C1*(1-nu),(C1*nu),0],[(C1*nu),C1*(1-nu),0],[0,0,(((1-(2*nu))/2)*C1)]]);     #Plane strain condition
  if Discontinuity_2 == 'Yes':
        C2=(E2/1-(nu2**2));
        D2=np.array([[C2,(C2*nu2),0],[(C2*nu2),C2,0],[0,0,(((1-nu2)/2)*C2)]]);   #Material Tangent matrix for an inclusion
  else:
        D2=np.array([[0,0,0],[0,0,0],[0,0,0]])

elif condition == 2:      ###Plane strain condition for 2
  C=E/((1-nu)*(1-(2*nu))) 
  Dstd=np.array([[C*(1-nu),(C*nu),0],[(C*nu),C*(1-nu),0],[0,0,(((1-(2*nu))/2)*C)]]);

  C1=E1/((1-nu1)*(1-(2*nu1))) 
  D1=np.array([[C1*(1-nu1),(C1*nu1),0],[(C1*nu1),C1*(1-nu1),0],[0,0,(((1-(2*nu1))/2)*C1)]]); 
   
  if Discontinuity_2 == 'Yes':
        C2=E2/((1-nu2)*(1-(2*nu2))) 
        D2=np.array([[C2*(1-nu2),(C2*nu2),0],[(C2*nu2),C2*(1-nu2),0],[0,0,(((1-(2*nu2))/2)*C2)]]);    #Material Tangent matrix for an inclusion
  else:
        D2=np.array([[0,0,0],[0,0,0],[0,0,0]])
#Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);
#Denr=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);

#########################Stiffness matrix for individual elements (K_local)############################
Bmatrix=[];
totDOF=usize+aDOFsiz;
#K=np.zeros((totDOF,totDOF));      
for i in range(len(element)):        
      Kenr=0;
      Kstd=0;
      #sigmaenr=0;
      #sigmastd=0;
     # sigmaenrfl=0;
      if i in elementc:    ##To create enriched B-matrix(Updated one)
        Dvalues3=[];
        for k in range(len(parlist[part])):
            Dvalues1=[];
            for j in range(4):
                  #Bmatenr,detstd,D= my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc,phivalforpart,phivalues); ##Old(01-04-23)
                  Bmatenr,detstd,detenrpart,D,Dvalues2= my_funcenh(i,j,parlist,part,k,Dstd,D1,D2,elementcuthole,elementcutinc,Dvalues1,
                                                                    xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc,phivalforpart);  ##original
                  KK2=np.matmul(np.transpose(Bmatenr),D);
                  #sigmaenr1=np.matmul(D,Bmatenr)*detstd
                  K2=wi[j]*np.matmul(KK2,Bmatenr)*detenrpart*detstd;    ###Detstd gives determinant for the single element
                  #K2=wi[j]*np.matmul(KK2,Bmatenr)*detstd; 
                  Kenr=K2+Kenr;
                  #sigmaenr=sigmaenr1+sigmaenr;
                  K4=(Kenr);
            Dvalues3.append(Dvalues2)
           # sigmaenrfl=(sigmaenr)+sigmaenrfl;                            ###multiplication of D and Bmatenr for all single subgrids
        Bmatrix.append(Bmatenr);
        #sigma.append(sigmaenrfl);
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
             #  sigmastd=sigmastd1+sigmastd;
               K4=(Kstd);
            Bmatrix.append(Bmatstd); 
            #sigma.append(sigmastd); 
              
           
      Kloc.append((K4));
     #sigma.append(KK2);  ###Sigma calc is wrong(Detected :12-04-23)
      det=detstd;
      
##################################################################################################


##################Assembly of K-matrix(Stiffness matrix)#######################
K=np.zeros((totDOF,totDOF));   ##Initialising the size of the global stiffness matrix

#function call (assembly_K) which assembles all the local stiffness of the elements in the global stiffness matrix
 
K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)  
###################################################################################


         
flatteig,flatt1=np.linalg.eig(Kloc[23]);#To check the whether the stiffness matrix is singular or not(it is a singular matrix)
flatteig1,flatt2=np.linalg.eig(Kloc[26]);
sum1=np.sum(flatteig)
sum2=np.sum(flatteig1)
#Flatt=np.linalg.solve(K,F)

##########################To build the element body force matrix:-###################################
##Body forces acting on the plate will be same for all the elements because forces acts equally
##on all elements under gravity.

bodyforceelem=[];   ##List used to store body forces for an element individually 

##Function call(Totbodyforcefunc) calculates the body forces of all the elements individually
#and save it in the list [bodyforceelem]
totbodyforce=Totbodyforcefunc(i,Bmatstd,Bmatenr,element,elementnoc,wi,t,det,bodforce,bodyforceelem,
                                      N1,N2,N3,N4,c1,c2,phivalues)


##Assembly of F matrix(body force)
Fbod1=np.zeros(np.shape(flatten));
#for i in range(len(element)):
   # fbloc1=totbodyforce[i];
    #Fbod assembles the 
Fbod=Totbodyforceassm(i,globnumber,elementnoc,DOF,globnodesextra,Fbod1,element,totbodyforce)
          
################################Traction force##################################################
           
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
     F2=[];
     #belement will store the element numbers which are located on the boundary whose input has been given
     #by the user in the form of Edge number
     belement=list(np.concatenate(belement1).flat) 
     F3=Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,wi,belement,
                     N1,N2,N3,N4,c1,c2,phivalues,Edge,l,F2)   
    #for i in range(len(element)):
     # fbloc2=F3[i];
     Ftract=Tottractforceassm(i,globnumber,elementnoc,DOF,globnodesextra,Ftract1,element,F3)  
   
  ii = input("Want to give Traction force for another edge Type(Yes or No): ")   #Yes or no
  if ii == "No":
    break; 

######################################################################################################

#######################______________'F'_vector______________######################################
###Adding global body force vector and global traction force vector and store it in 'F' variable

F=Fbod+Ftract;

#####################################################################################################

####################_______Dirichlet_Boundary_Conditions_________###################################


######First finding the fixed nodes(global number) which helps us to solve the linear system of 
######equations

BCnodes=[]; ##Checking for boundary nodes
fixedDOFSlist1=[];
Constraints=input("Are you need to fix/constraints the Edges? Type:(Yes or No)")
if Constraints == 'No':
    fixedDOFSlist1=[]
else:
 while True: 
  fixed_edge=input("Enter the name of the fixed edge : ")
  boundarynodes=[]
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
  
###################Applying Dirchelet Boundary Conditions###################
###After finding the fixed nodes (by knowing which edge is fixed),we can able to transform the Stiffness matrix(K)
###in order to solve the linear system of equations

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
        K,l3,F=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3,F)
        
  ii = input("Want to fix another edge,Please Type Yes or No : ")   #Yes or no
  if ii == "No":
      break;

####Temp#####################################################################
fixedDOFSlist=['x'+str(0),'y'+str(0),'x'+str(1640),'y'+str(1640)]#,\]
               #'a'+str(832),'b'+str(832),'x'+str(848),'y'+str(848),'a'+str(848),'b'+str(848)]
                  # 'x'+str(1168),'y'+str(1168),'a'+str(1168),'b'+str(1168)];
for i in range(len(fixedDOFSlist)):
      fixDOF=fixedDOFSlist[i];
      K,l3,F=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3,F)
##############################################################################
      

##############___________Solving tne Linear system of equations________________############################
###Using linalg from numpy library,solve the linear system of equation to get nodal displacements
###of all the nodes.

flatten=np.linalg.solve(K,F);  

#######################################################################################################

######################_________Results______________#######################################################
###To calculate displacements,stresses,strains for all elements:-
increment=0;
for i in range(len(element)):
    disp=umatrix[i];   
    disp3,actdisp=Disp_stress_strain(i,disp,flatten,Bmatrix,elementc,elementnoc,elementinsidehole,phivalues)         ##Stress and strain will come from this block
    displacement.append(actdisp);
    stress1,strain1=stress_strain_component(i,disp,flatten,element,elementc,elementnoc,elementinsidehole,Dstd,D1,phivalues,\
                                            phivalforpart,increment)
    stress.append(stress1)
    strain.append(strain1)
 

 
xvalues2=np.linspace(0,br,nor)   ###Arranging elements
yvalues2=np.linspace(0,L,noc)    ###Arranging elements
xv1, yv1 = np.meshgrid(xvalues2, yvalues2)
zvalues2=np.zeros([noc,nor])
zvalues3=np.zeros([noc,nor])
loop=0;
for i in range(len(yvalues2)):
    for j in range(len(xvalues2)):
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
zvalues=np.zeros([noc,nor])
zvalues1=np.zeros([noc,nor])
loop1=0;
for i in range(len(yvalues)):
    for j in range(len(xvalues)):
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
plt.savefig('Books_read.png')  ##Working
plt.show()
plt.show()       
