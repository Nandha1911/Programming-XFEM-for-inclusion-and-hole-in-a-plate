import numpy as np
import math
import pytest
import sympy as sp
from Building_elements import creatingelement
from LevelSetFunction import LevelSet
from Rectangular_sub_grid_nodes import Partitioning_method
from Determinant_function import my_funcdet

#Gauss points
xi1 = np.array([(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3),math.sqrt(1/3)]);
xi2 = np.array([(math.sqrt(1/3)),(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3)]); 
#Shape Functions:-
c1 = sp.Symbol('c1');
c2 = sp.Symbol('c2');

N1=((1-c1)*(1+c2))/4;
N2=((1-c1)*(1-c2))/4;      
N3=((1+c1)*(1-c2))/4;
N4=((1+c1)*(1+c2))/4;
N=np.array([N1,N2,N3,N4]);


 
#class test_x(unittest.TestCase):
def test_coordinates():
        i=1;j=0;l=0.1
        result=np.array(creatingelement(i,j,l))
        expected_output=np.array([(0,0.1),(0,0.2),(0.1,0.2),(0.1,0.1)])
        print ("\n")
        for i in range(len(result)):
            if (result[i] == expected_output[i]).all():
                  print ("The coordinates are correct")
                  print ("The results are: ",result[i])
                  print ("The expected output: ",expected_output[i])
                  print ("\n")
                  assert True
            else:
                  print ("The coordinates are incorrect")
                  print ("\n")
                  assert False
     
def test_curveintersect(): #498,1269,539,1748,149
                a1=2;b1=2;r1=0.8;a2=3;b2=4;r2=0.8;
                phole=[];pinc=[];
                testelem=np.array([[(1.8,1.2),(1.8,1.3),(1.9,1.3),(1.9,1.2)],
                       [(2.9,3.1),(2.9,3.2),(3,3.2),(3,3.1)],[(1.9,1.3),(1.9,1.4),(2.0,1.4),(2.0,1.3)], 
                        [(2.8,4.3),(2.8,4.4),(2.9,4.4),(2.9,4.3)],[(3.6,0.2),(3.6,0.3),(3.7,0.3),(3.7,0.2)]])
                exp_value=np.array([[1,-1,-1,1],[1,1,-1,1],[-1,-1,-1,-1],[-1,-1,-1,-1],[1,1,1,1]])
                print ("\n")
                for i in range(len(testelem)):
                    phole=LevelSet(testelem[i],a1,b1,r1)
                    pinc=LevelSet(testelem[i],a2,b2,r2)
                    case='case'+str(i)
                
                    if case == 'case0':
                        if np.array_equal(phole,exp_value[0])==True and np.array_equal(pinc,exp_value[0])==False:
                           print ("The curve(hole) passes through this element(498)")
                           print ("Output value",phole)
                           print ("Expected value",exp_value[0])
                           print ("The curve(inclusion) does not pass through this element(498)")
                           print ("\n")
                           assert True
                        else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False
                  
                    elif case == 'case1':
                        if np.array_equal(phole,exp_value[1])==False and np.array_equal(pinc,exp_value[1])==True:
                            print ("The curve(inclusion) pass through this element(1269)")
                            print ("Output value",pinc)
                            print ("Expected value",exp_value[1])
                            print ("The curve(hole) does not pass through this element(1269)")
                            print ("\n")
                            assert True 
                        else:
                            print ("The level set function is wrongly calculated for this coordinates")
                            assert False

                    elif case == 'case2':
                        if np.array_equal(phole,exp_value[2])==True and np.array_equal(pinc,exp_value[2])==False:
                           print ("The element(539) is inside the hole")
                           print ("Output value",phole)
                           print ("Expected value",exp_value[2])
                           print ("The element(539) is not inside the inclusion")
                           print ("\n")
                           assert True
                        else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False

                    elif case == 'case3':
                        if np.array_equal(phole,exp_value[3])==False and np.array_equal(pinc,exp_value[3])==True:
                           print ("The element(1748) is inside the inclusion")
                           print ("Output value",pinc)
                           print ("Expected value",exp_value[3])
                           print ("The element(1748) is not inside the hole")
                           print ("\n")
                           assert True
                        else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False

                    elif case == 'case4':
                         if np.array_equal(phole,exp_value[4])==True and np.array_equal(pinc,exp_value[4])==True:
                             print ("Both the curve(hole and inclusion)does not pass through this element(149)")
                             print ("Output value",phole,pinc)
                             print ("Expected value",exp_value[4])
                             print ("\n")
                             assert True
                         else:
                             print ("The level set function is wrongly calculated for this coordinates")
                             assert False
                       

def test_partitionelement():
             i=0;l=0.1;elementcc=i;element1=XFEM_python_32.element #element=np.array([(1.9,1.1),(1.9,1.2),(2,1.2),(2,1.1)])
             exp_value=np.array([[(1.9,1.1),(1.9,1.13),(1.93,1.13),(1.93,1.1)],[(1.93,1.1),(1.93,1.13),(1.96,1.13),(1.96,1.1)],[(1.96,1.1),(1.96,1.13),(2,1.13),(2,1.1)],
                                 [(1.9,1.13),(1.9,1.16),(1.93,1.16),(1.93,1.13)],[(1.93,1.13),(1.93,1.16),(1.96,1.16),(1.96,1.13)],[(1.96,1.13),(1.96,1.16),(2,1.16),(2,1.13)],
                                 [(1.9,1.16),(1.9,1.2),(1.93,1.2),(1.93,1.16)],[(1.93,1.16),(1.93,1.2),(1.96,1.2),(1.96,1.16)],[(1.96,1.16),(1.96,2),(2,1.2),(2,1.16)]])
             output_value=Partitioning_method(i,l,elementcc,element1)
             print ("\n")
             print (output_value)
             if np.array_equal(exp_value,output_value):
                     print("The expected and output values are same")
                     assert True
             else:
                     print("The partitioning function is calculating wrongly")
                     assert False 
#if __name__ =='__main__':
#    unittest.main()
 

##################Sanity Test###########################
def test_shapefunctions():
    ##Add the shape functions N1,N2,N3,N4 and the expected result should be 1 which proves that the substituted shape functions in the XFEM program is  correct
    expected_val=1;
    #Gauss points(xi1,xi2)
    # xi1 = np.array([(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3),math.sqrt(1/3)]);
     #xi2 = np.array([(math.sqrt(1/3)),(-math.sqrt(1/3)),(-math.sqrt(1/3)),math.sqrt(1/3)]); 
    #Shape Functions:-
     #c1 = sp.Symbol('c1');
     #c2 = sp.Symbol('c2');

    #N1=((1-c1)*(1+c2))/4;
    #N2=((1-c1)*(1-c2))/4;      
    #N3=((1+c1)*(1-c2))/4;
    #N4=((1+c1)*(1+c2))/4;
    #N=np.array([N1,N2,N3,N4]);
    print ("\n")
    print ("To check the sum of shape functions with respect to their Gauss points are one")
    print ("\n")
    for j in range(len(N)):
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
    for j in range(len(N)):
        dN1dxi1=N1.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);   #Differentation w.r.t c1
        dN2dxi1=N2.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
        dN3dxi1=N3.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
        dN4dxi1=N4.diff(c1).subs(c1,xi1[j]).subs(c2,xi2[j]);
        dN1dxi2=N1.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);   #Differentation w.r.t c2
        dN2dxi2=N2.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
        dN3dxi2=N3.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
        dN4dxi2=N4.diff(c2).subs(c1,xi1[j]).subs(c2,xi2[j]);
        Gauss_point=1+j;
        output_val=dN1dxi1+dN2dxi1+dN3dxi1+dN4dxi1+dN1dxi2+dN2dxi2+dN3dxi2+dN4dxi2
        if expected_val == output_val:
          print ('Output value is',output_val)
          print ('The sum of derivative of shape functions with respect to Gauss point',Gauss_point,'is equal to zero')
          print ("\n")
          assert True
        else:
           print ('The sum of derivative of shape functions with respect to Gauss point',Gauss_point,'is not  equal to zero')
           print ("\n")
           assert False
        
def test_Jacobimatrix():
      print ("\n")
      for j in range(4):      ##This j loop will run for 4 Gauss points
    
           i=0;element=[np.array([[0,0],[0,1],[1,1],[1,0]])];elementnoc=[0]
           dNdxi,detj,Jac=my_funcdet(i,j,xi1,xi2,c1,c2,N1,N2,N3,N4,element,elementnoc)
           Theta=115         ###Theta value given for the rotation of Jacobi matrix
           expected_val=detj
           print ('Jacobi matrix before rotation')
           print (Jac)
           print ('Expected_val =',expected_val)  ####Determinant of Jacobi
           Q=np.array([[math.cos(Theta),math.sin(Theta)],[-math.sin(Theta),math.cos(Theta)]])
           Rot1=np.matmul(np.transpose(Q),Jac)
           RotJacobi=np.float32(np.matmul(Rot1,Q))    ###Determinant of Jacobi after Rotation
           print ('Jacobi matrix after rotation')
           print (RotJacobi)
           if np.linalg.det(RotJacobi) == expected_val:
                 output_val=np.linalg.det(RotJacobi);
                 print ("Output value=",output_val)
                 print ("The determinant of jacobi matrix after rotation is equal to the determinant of Jacobi matrix before rotation with respect to Gauss point",j+1)
                 print ("\n")
                 assert True
           else:
                 output_val=np.linalg.det(RotJacobi);
                 print ("Output value=",output_val)
                 print ("The determinant of jacobi matrix after rotation is not equal to the determinant of Jacobi matrix before rotation with respect to Gauss point",j+1)
                 print ("\n")
                 assert False
