import numpy as np
import math
import sympy as sp
#from Building_elements import creatingelement 
from Naming_the_umatrix import naming_umatrix
from B_matrix_for_std import my_func
from Assembly_Stiffness_mat import assembly_K
from Building_elements_1 import creatingelement
from Dirichlet_Boundarycondition import DirichletBC
from Boundary_elements import Tractionforce
from Traction_force_function_1 import Tractforce
from Traction_force_assembly import Tottractforceassm

#Input:-
L=2;
b=2;
l=1;nor= int((b/l)); noc= int((L/l));
element=[];
globnodes=[];
globU=[];
u=[];
aDOF=[];
globnodesextra=[];
globnumber=[];
Kloc=[];
umatrix=[];
elementswglobnode=[];
fixedDOFSlist1=[];
l3=[]; ##Checking
belement1=[];
DOF=2
Bmatrix=[];
t=1;
F3=[];

#Material_properties:-
E1=3*(10**7) 
E2=0.01*E1;            #Young's modulus in N/cm2(Hole) #To avoid singular stiffness matrices1(E2/E1=0.01)
E3=4.19*(10**7)
nu=0.3

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

for i in range(noc):  
    for j in range(nor):
        element.append(creatingelement(i,j,l))

#element1=np.array([[0,1],[0,0],[2,0.5],[2,1]])
#element.append(element1)    
elementnoc=[0,1,2,3]
elementc=[];         ##Elements that are cut by the curve has not been considered
elementinsidehole=[]
elementinsideinc=[]
new=[];

zz=0;
aDOFsiz=0;
usize=0;

globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF=naming_umatrix(noc,nor,l,new,zz,aDOFsiz,usize,globnodesextra,globnodes,globU,u,aDOF)

flatten1 = list(np.concatenate(globU).flat); 
flatten = np.array([flatten1]).T;    

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

#globnumber1=np.array([0,1,2,3])  
#globnumber.append(globnumber1) 
    
C=(E1/1-(nu**2));
Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate

C1=(E2/1-(nu**2));
D1=np.array([[C1,(C1*nu),0],[(C1*nu),C1,0],[0,0,(((1-nu)/2)*C1)]]);  #Material Tangent Matrix for a hole

C2=(E3/1-(nu**2));
D2=np.array([[C2,(C2*nu),0],[(C2*nu),C2,0],[0,0,(((1-nu)/2)*C2)]]);   #Material Tangent matrix for an inclusion
totDOF=usize+aDOFsiz;
for i in range(len(element)): 
    Kstd=0;
    for j in range(4): ##To create standard B-matrix
       #J,Bmat= my_func(i,j);
       #Jac=np.matmul(J,element[i]);
       Bmatstd,detstd= my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
       if i in elementinsidehole:##New line
           KK2=np.matmul(np.transpose(Bmatstd),D1);
       elif i in elementinsideinc:
           KK2=np.matmul(np.transpose(Bmatstd),D2);
       else:   
           KK2=np.matmul(np.transpose(Bmatstd),Dstd);
       K1=wi[j]*np.matmul(KK2,Bmatstd)*detstd;
       Kstd=K1+Kstd;
       K4=(Kstd);
    Bmatrix.append(Bmatstd)
       
    Kloc.append(K4)
    
K=np.zeros((totDOF,totDOF));
K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)
eigenvalue1,eigenvector=np.linalg.eig(K)
eigenvaluewbound=np.round(eigenvalue1,4)
unconstrained_rigid_body_motions=len(eigenvaluewbound[eigenvaluewbound==0])
print (unconstrained_rigid_body_motions)

##This block(166 to 199) has been written in order to find the constrained_rigid_body_motions(Total Eigen value=0)
'''
FixedDOFS3=input("Do you need to constrain the nodes of respective edge in the x-dir or y-dir or both?")
#boundarynodes=[20,15,10,5,0]
boundarynodes=[6,3,0]
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
   Kbound,l3=DirichletBC(i,fixDOF,fixedDOFSlist,K,flatten1,l3)
   
eigenvalue2,eigenvector = np.linalg.eig(Kbound);
eigenvaluebound=np.round(eigenvalue2,4)
constrained_rigid_body_motions=len(eigenvaluebound[eigenvaluebound==0])
print (constrained_rigid_body_motions)
'''
'''
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
belement=list(np.concatenate(belement1).flat)

Bmatenr=np.array([0]);
elementc=[];phivalues=[];
for i in range(len(element)):    ###Error in Traction_force_function_1(Last Update:14/04/23)
   F2=Tractforce(i,Bmatstd,Bmatenr,element,elementc,elementnoc,tractforces,t,wi,belement,
                 N1,N2,N3,N4,c1,c2,xi1,xi2,phivalues,Edge)
   F3.append(F2);

Ftract1=np.zeros(np.shape(flatten));
for i in range(len(element)):
  fbloc2=F3[i];
  Ftract=Tottractforceassm(i,globnumber,elementnoc,DOF,globnodesextra,fbloc2,Ftract1)  
'''
'''
####This block is return for calculating Rigid body displacement test(234 to 257)
flatten1=np.array([2,2,2,2,2,2,2,2,0,0,2,2,2,2,2,2,2,2])
flatten = np.array([flatten1]).T; 
#Rf=np.matmul(K,flatten)
Rf=K@flatten
F1=np.zeros([18,1])
F=F1-Rf
for i in range(18):
    for j in range(18):
        if i == 8 and j==8 or i==9 and j==9 or i==8 and j==9 or i==9 and j==8:
            break;
        elif i == j:
            K[i,j]=1;
        else:
            K[i,j]=0;
            K[j,i]=0;
            
for i in range(18):
    if i==8 or i== 9:
        F[i]=F[i]
    else:
        F[i]=0;
        
            
U=np.linalg.solve(K,F)
####This block is return for calculating Rigid body displacement test(234 to 257)
'''
Theta=45;
Q=np.array([[math.cos(Theta),math.sin(Theta)],[-math.sin(Theta),math.cos(Theta)]])
flatten1=np.zeros([18,1])

for i in range(len(globnodes)):
    if i==4:
        rotatedpos=np.matmul(np.array([0,0]),Q)
        flatten1[(i*2),0]=rotatedpos[0]
        flatten1[((i*2)+1),0]=rotatedpos[1]
    
    else:
        rotatedpos=np.matmul(globnodes[i],Q)
        flatten1[(i*2),0]=rotatedpos[0]
        flatten1[((i*2)+1),0]=rotatedpos[1]
        
flatten = flatten1
Rf=K@flatten
F1=np.zeros([18,1])
F=F1-Rf 

print(globnodes[4])
newpoints=np.matmul(globnodes[4],Q)

for i in range(18):
    for j in range(18):
        if i == 8 and j==8 or i==9 and j==9 or i==8 and j==9 or i==9 and j==8:
            break;
        elif i == j:
            K[i,j]=1;
        else:
            K[i,j]=0;
            K[j,i]=0;
            
for i in range(18):
    if i==8 or i== 9:
        F[i]=F[i]
    else:
        F[i]=0;

U=np.linalg.solve(K,F)
tol=1e-5
print (abs(U[8]-newpoints[0]))
print (abs(U[9]-newpoints[1]))
if abs(U[8]-newpoints[0]) < tol and abs(U[9]-newpoints[1]) < tol:
    print ("Test case passed")