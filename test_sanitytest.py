import numpy as np
import math
import sympy as sp
import pytest
 
from Naming_the_umatrix import naming_umatrix
from Determinant_function import my_funcdet
from B_matrix_for_std import my_func
from Assembly_Stiffness_mat import assembly_K
from Building_elements_1 import creatingelement
from Dirichlet_Boundarycondition import DirichletBC
from Boundary_elements import Tractionforce
from Traction_force_function_2 import Tractforce
from Traction_force_assembly import Tottractforceassm
from glob_number_allocation import glob_number
from Boundary_nodes import fixednodes

#Input:-
L=2;
b=2;
l=1;nor= int((b/l)); noc= int((L/l));
#element=[];
#DOF=2
#Bmatrix=[];


#Material_properties:-
E=1*(10**5) 
E1=0.01*E;            #Young's modulus in N/cm2(Hole) #To avoid singular stiffness matrices1(E2/E1=0.01)
E2=4.19*(10**5)
nu=0.3;nu1=0.3;nu2=0.3;
##Gauss points:-
xi1 = np.array([(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3),math.sqrt(1/3)]);
xi2 = np.array([(-math.sqrt(1/3)),(math.sqrt(1/3)),(-math.sqrt(1/3)),(math.sqrt(1/3))]);
##Weight_functions
wi = np.array([1,1,1,1]);


c1 = sp.Symbol('c1');
c2 = sp.Symbol('c2');
#Shape Functions:-
N1=((1-c1)*(1-c2))/4;
N2=((1+c1)*(1-c2))/4;      
N3=((1+c1)*(1+c2))/4;
N4=((1-c1)*(1+c2))/4;

Nstd=np.array([[N1,0,N2,0,N3,0,N4,0],[0,N1,0,N2,0,N3,0,N4]])


##################Sanity Test###########################
def test_shapefunctions():
    ##Add the shape functions N1,N2,N3,N4 and the expected result should be 1 which proves that the substituted shape functions in the XFEM program is  correct
    expected_val=1;
    #Gauss points(xi1,xi2)
    print ("\n")
    print ("To check the sum of shape functions with respect to their Gauss points are one")
    print ("\n")
    for j in range(4):  ## 4 Gauss points ##
       N11=N1.subs(c1,xi1[j]).subs(c2,xi2[j])
       N12=N2.subs(c1,xi1[j]).subs(c2,xi2[j])
       N13=N3.subs(c1,xi1[j]).subs(c2,xi2[j])
       N14=N4.subs(c1,xi1[j]).subs(c2,xi2[j])
       Gauss_point=1+j;
       output_val=N11+N12+N13+N14
       if expected_val == output_val:
          print ('Output value is',output_val)
          print ('The sum of shape functions with respect to Gauss point',Gauss_point,'is equal to one')
          print ("\n")
          assert True
       else:
          print ('The sum of shape functions with respect to Gauss point',Gauss_point,'is not  equal to one')
          print ("\n")
          assert False
      


def test_derivativeshapefun():
    expected_val=0;
    print ("\n")
    print ("To check the sum of derivative of shape functions with their respective Gauss points are zero")
    print ("\n")
    for j in range(4):  ##4 Gauss points##
        dNdxi1=np.array([[(xi2[j]-1)/4,(1-xi2[j])/4,(1+xi2[j])/4,(-xi2[j]-1)/4],[(xi1[j]-1)/4,(-xi1[j]-1)/4,(1+xi1[j])/4,(1-xi1[j])/4]])
        Gauss_point=1+j;
        output_val=np.sum(dNdxi1)
        if expected_val == output_val:
          print ('Output value is',output_val)
          print ('The sum of derivative of shape functions with respect to Gauss point',Gauss_point,'is equal to zero')
          print ("\n")
          assert True
        else:
           print ('The sum of derivative of shape functions with respect to Gauss point',Gauss_point,'is not  equal to zero')
           print ("\n")
           assert False
           
#####################        
def test_Jacobimatrix():
         print ("\n")
         for j in range(4):      ##This j loop will run for 4 Gauss points
           i=0;element=[np.array([[0,1],[0,0],[1,0],[1,1]])];elementnoc=[0]
           dNdxi,detj,Jac=my_funcdet(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc)
           Theta=115         ###Theta value given for the rotation of Jacobi matrix
           expected_val=detj
           print ('Jacobi matrix before rotation')
           print ('Expected_val =',expected_val)  ####Determinant of Jacobi
           Q=np.array([[math.cos(Theta),math.sin(Theta)],[-math.sin(Theta),math.cos(Theta)]])
           Rot1=np.matmul(np.transpose(Q),Jac)
           Rotjacobi=round(np.linalg.det(Rot1),4)
           print ('Jacobi matrix after rotation')
           print (Rotjacobi)
           if  Rotjacobi == expected_val:
                 print ("Output value=",Rotjacobi)
                 print ("The determinant of jacobi matrix after rotation is equal to the determinant of Jacobi matrix before rotation with respect to Gauss point",j+1)
                 print ("\n")
                 assert True
           else:
                 print ("Output value=",Rotjacobi)
                 print ("The determinant of jacobi matrix after rotation is not equal to the determinant of Jacobi matrix before rotation with respect to Gauss point",j+1)
                 print ("\n")
                 assert False
                 

def test_eigenvaluetest():
  print("Sanity_test_4a")
  element=[]
  expectedeigvalforunconstrained=3;expectedeigvalforconstrained=0;                
  for i in range(noc):  
    for j in range(nor):
        elementcall=creatingelement(i,j,l)
        element.append(elementcall)

    
  elementnoc=[0,1,2,3]
  elementc=[];         ##Elements that are cut by the curve has not been considered
  new=[];

  zz=0;
  aDOFsiz=0;
  usize=0;globnodesextra=[];globU=[];umatrix=[];l3=[];elementswglobnode=[];Kloc=[];l3=[];DOF=2;
  Bmatrix=[]

  globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF,cornernodenum,allboundarynodenum,ordinarynodenum=naming_umatrix(noc,nor,l,new,aDOFsiz,\
                                                                                                          usize,globnodesextra,globU)

  flatten1 = list(np.concatenate(globU).flat); 
  flatten = np.array([flatten1]).T;    
      
  count=0;  
  elementswglobnode,globnumber=glob_number(noc,nor,count,elementc,u,aDOF) 
    
  C=(E/1-(nu**2));
  Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate

  totDOF=usize+aDOFsiz;
  for i in range(len(element)): 
    Kstd=0;
    for j in range(4): ##To create standard B-matrix
       Bmatstd,detstd= my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
       KK2=np.matmul(np.transpose(Bmatstd),Dstd);
       K1=wi[j]*np.matmul(KK2,Bmatstd)*detstd;
       Kstd=K1+Kstd;
       K4=(Kstd);
    Bmatrix.append(Bmatstd)
       
    Kloc.append(K4)
    
  K=np.zeros((totDOF,totDOF));
  K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)
  eigenvalue1,eigenvector=np.linalg.eig(K)
  eigenvalueunconstrained=np.round(eigenvalue1,4)
  unconstrained_rigid_body_motions=len(eigenvalueunconstrained[eigenvalueunconstrained==0])
  print ('unconstrained_rigid_body_motions = ',unconstrained_rigid_body_motions)
  if unconstrained_rigid_body_motions ==  expectedeigvalforunconstrained:
      print("3 zero eigen values for unconstrained rigid body motions.So the nature of Stiffness matrix is conserved")
      assert True
  else:
      print("The nature of Stiffness matrix is violated")
      assert False

  fixedDOFSlist1=[];F=np.zeros(np.shape(flatten));
  print("Sanity_test_4b")
  print('######################################################')  
  print("Applying Dirichlet Boundary Conditions")
  fix="1"
  boundarynodes=[]
  if fix == '1':
    fixed_edge='Edge1'   ###Edge1 will be fixed
    #boundarynodes=[]
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
   ##############Cornernodes##########################################
  elif fix == '2':
     corner_node=(input("Enter the number of the cornernode need to be fixed (Type 1 or 2 or 3 or 4) : "))
     if corner_node == '1':   
          bnodes = [cornernodenum[0]];
     elif corner_node == '2':   
          bnodes = [cornernodenum[1]]; 
     elif corner_node == '3':   
         bnodes = [cornernodenum[2]];
     elif corner_node == '4':   
          bnodes = [cornernodenum[3]]; 
     else:
       print("The entered edge has not been found amoung those edges")
     
  boundarynodes.extend(bnodes)   ##Bnodes and boundary nodes are same
  
###################Applying Dirichlet Boundary Conditions###################
###After finding the fixed nodes (by knowing which edge is fixed),we can able to transform the Stiffness matrix(K)
###in order to solve the linear system of equations

  FixedDOFS3="both"     ###Both x and y -axis has been fixed
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

  eigenvalue2,eigenvector = np.linalg.eig(K);
  eigenvalueconstrained=np.round(eigenvalue2,4)
  constrained_rigid_body_motions=len(eigenvalueconstrained[eigenvalueconstrained==0])
  print ('constrained_rigid_body_motions = ',constrained_rigid_body_motions)
  if constrained_rigid_body_motions == expectedeigvalforconstrained:
      print("zero eigen values for constrained rigid body motions.So the nature of Stiffness matrix is conserved")
      assert True
  else:
      print("The nature of Stiffness matrix is violated")
      assert False


####This block is return for calculating Rigid body displacement test(234 to 257)
def test_rigidbodydisplacementtest():
  print("Sanity_test_5")
  element=[]
  for i in range(noc):  
    for j in range(nor):
        elementcall=creatingelement(i,j,l)
        element.append(elementcall)
    
  elementnoc=[0,1,2,3]
  elementc=[];         ##Elements that are cut by the curve has not been considered
  #elementinsidedis1=[]
  #elementinsidedis2=[]
  new=[];

  zz=0;
  aDOFsiz=0;
  usize=0;globnodesextra=[];globU=[];umatrix=[];elementswglobnode=[];Kloc=[];DOF=2;
  Bmatrix=[];

  globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF,cornernodenum,allboundarynodenum,ordinarynodenum=naming_umatrix(noc,nor,l,new,aDOFsiz,\
                                                                                                          usize,globnodesextra,globU)
   
  count=0;  
  elementswglobnode,globnumber=glob_number(noc,nor,count,elementc,u,aDOF) 
  C=(E/1-(nu**2));
  Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate
 
  expected_result=[2,2]
  #####Adding displacement=2 to all the boundary nodes except center node#####
  flatten1=np.array([2,2,2,2,2,2,2,2,0,0,2,2,2,2,2,2,2,2])  
  flatten = np.array([flatten1]).T; 

  #####Now calculating individual stiffness matrices##########
  totDOF=usize+aDOFsiz;
  for i in range(len(element)): 
   Kstd=0;
   for j in range(4): ##To create standard B-matrix
     #J,Bmat= my_func(i,j);
     #Jac=np.matmul(J,element[i]);
     Bmatstd,detstd= my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
     KK2=np.matmul(np.transpose(Bmatstd),Dstd);
     K1=wi[j]*np.matmul(KK2,Bmatstd)*detstd;
     Kstd=K1+Kstd;
     K4=(Kstd);
  Bmatrix.append(Bmatstd)
     
  Kloc.append(K4)
  
  K=np.zeros((totDOF,totDOF));
  K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)
  Residualforces=np.matmul(K,flatten)
  F1=np.zeros([18,1])
  F=F1-Residualforces
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
  if U[8] == expected_result[0] and U[9] == expected_result[1]:
    print("Rigid body displacement test has been successfully done.The centre node will also move according to the given displacement to all the boundary nodes")
    assert True
  else:
   print("Rigid body displacement test has been failed")
   assert False

####This block is return for calculating Rigid body displacement test(234 to 257)
def test_rigidbodyrotationtest():
  print("Sanity_test_6")
  element=[]
  for i in range(noc):  
    for j in range(nor):
      elementcall=creatingelement(i,j,l)
      element.append(elementcall)
  
  elementnoc=[0,1,2,3]
  elementc=[];         ##Elements that are cut by the curve has not been considered
  new=[];

  zz=0;
  aDOFsiz=0;
  usize=0;globnodesextra=[];globU=[];umatrix=[];elementswglobnode=[];Kloc=[];DOF=2;
  Bmatrix=[]

  globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF,cornernodenum,allboundarynodenum,ordinarynodenum=naming_umatrix(noc,nor,l,new,aDOFsiz,\
                                                                                                        usize,globnodesextra,globU)
 
  count=0;  
  elementswglobnode,globnumber=glob_number(noc,nor,count,elementc,u,aDOF) 
  C=(E/1-(nu**2));
  Dstd=np.array([[C,(C*nu),0],[(C*nu),C,0],[0,0,(((1-nu)/2)*C)]]);  ##Material Tangent matrix for plate


#####Now calculating individual stiffness matrices##########
  totDOF=usize+aDOFsiz;
  for i in range(len(element)): 
    Kstd=0;
    for j in range(4): ##To create standard B-matrix
      Bmatstd,detstd= my_func(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc);
      KK2=np.matmul(np.transpose(Bmatstd),Dstd);
      K1=wi[j]*np.matmul(KK2,Bmatstd)*detstd;
      Kstd=K1+Kstd;
      K4=(Kstd);
  Bmatrix.append(Bmatstd)
   
  Kloc.append(K4)

  K=np.zeros((totDOF,totDOF));
  K,umatrix=assembly_K(K,Kloc,element,elementnoc,globnumber,DOF,globnodesextra,umatrix)
  Theta=135;
  Q=np.array([[math.cos(Theta),math.sin(Theta)],[-math.sin(Theta),math.cos(Theta)]])
  flatten=np.zeros([18,1])

  for i in range(len(globnodes)):
    if i==4:   ####Except centre node 4,all other nodes will be multiplied with rotation matrix
        rotatedpos=np.matmul(np.array([0,0]),Q)
        flatten[(i*2),0]=rotatedpos[0]
        flatten[((i*2)+1),0]=rotatedpos[1]
    
    else:
        rotatedpos=np.matmul(globnodes[i],Q)
        flatten[(i*2),0]=rotatedpos[0]
        flatten[((i*2)+1),0]=rotatedpos[1]
        
  Residualforces=np.matmul(K,flatten)
  F1=np.zeros([18,1])
  F=F1-Residualforces 

  print('centre_node =',globnodes[4])
  rotatedcentrenode=np.matmul(globnodes[4],Q)
  print('After centre node rotation =',rotatedcentrenode)

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

  if (U[8]-rotatedcentrenode[0]) < tol and (U[9]-rotatedcentrenode[1]) < tol:
    print("Rigid body rotation test has been successfully done.The centre node will also rotate according to the given displacement to all the boundary nodes")
    assert True
  else:
    print("Rigid body rotation test has been failed.")
    assert False
