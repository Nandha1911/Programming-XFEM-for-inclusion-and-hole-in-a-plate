import numpy as np
import math
import pytest
import sympy as sp
from Building_elements_1 import creatingelement
from LevelSetFunction import LevelSet
from Rectangular_sub_grid_nodes_1 import Partitioning_method
from Naming_the_umatrix import naming_umatrix
from glob_number_allocation import glob_number


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

N=np.array([N1,N2,N3,N4]);



##################For testing whether the coordinates of an element are generated correctly##################
def test_coordinates():
        print ("Unit_Test_1")
        i=0;j=1;l=1
        result=np.array(creatingelement(i,j,l))
        expected_output=np.array([(1,1),(1,0),(2,0),(2,1)])
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

##################For testing whether the element is standard or enriched element##################     
def test_levelsetmethod():
##Taking some random elements from 10*10 metre plate with mesh size of  1 metre

  print ("Unit_Test_2")
  a1=5;b1=5;r1=2.5;   ###a1 and b1 are the coordinates of the discontinuity where r1 is the radius of the
                    ### discontinuity
  discontinuity=[];
  testelem=np.array([[(0,3),(0,2),(1,2),(1,3)],
        [(3,3),(3,2),(4,2),(4,3)],[(4,3),(4,2),(5,2),(5,3)], 
         [(4,4),(4,3),(5,3),(5,4)]])  ###elements
  exp_value=np.array([[1,1,1,1],[1,1,1,-1],[-1,1,1,-1],[-1,-1,-1,-1]])
  print ("\n")

  for i in range(len(testelem)):
    discontinuity=LevelSet(testelem[i],a1,b1,r1)
    case='case'+str(i)
                
    if case == 'case0':
         if np.array_equal(discontinuity,exp_value[0])==True :
                           print ("The curve(hole) does not passes through this element")
                           print ("Output value",discontinuity)
                           print ("Expected value",exp_value[0])
                           print ("\n")
                           assert True
         else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False
                  
    elif case =='case1':
        if np.array_equal(discontinuity,exp_value[1])==True :
                            print ("The curve pass through this element")
                            print ("Output value",discontinuity)
                            print ("Expected value",exp_value[1])
                            print ("\n")
                            assert True 
        else:
                            print ("The level set function is wrongly calculated for this coordinates")
                            assert False

    elif case =='case2':
        if np.array_equal(discontinuity,exp_value[2])==True:
                           print ("The curve pass through this element")
                           print ("Output value",discontinuity)
                           print ("Expected value",exp_value[2])
                           print ("\n")
                           assert True
                           
        else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False

    elif case =='case3':
       if np.array_equal(discontinuity,exp_value[3])==True:
                           print ("The element is inside the curve")
                           print ("Output value",discontinuity)
                           print ("Expected value",exp_value[3])
                           print ("\n")
                           assert True
                         
       else:
                           print ("The level set function is wrongly calculated for this coordinates")
                           assert False


####################################For testing the rectangular-subgrids methods####################################
def test_partitionelement():
     print ("Unit_Test_3")
     i=0;l=1;elementc=[0];element=[np.array([[3,3],[3,2],[4,2],[4,3]])] 
     exp_value=[[[3,2.333],[3,2],[3.333,2],[3.333,2.333]],[[3.333,2.333],[3.333,2],[3.666,2],[3.666,2.333]],[[3.666,2.333],[3.666,2],[4,2],[4,2.333]],
                                 [[3,2.666],[3,2.333],[3.333,2.333],[3.333,2.666]],[[3.333,2.666],[3.333,2.333],[3.666,2.333],[3.666,2.666]],[[3.666,2.666],[3.666,2.333],[4,2.333],[4,2.666]],
                                 [[3,3],[3,2.666],[3.333,2.666],[3.333,3]],[[3.333,3],[3.333,2.666],[3.666,2.666],[3.666,3]],[[3.666,3],[3.666,2.666],[4,2.666],[4,3]]]
     output_value=Partitioning_method(i,l,elementc,element) 
     print ("\n")
     if np.allclose(exp_value,output_value):
                     print ("\n")
                     print("The expected subgrid nodes and output subgrid nodes are same")
                     assert True
                     print ("\n") 
     else:
                     print("The partitioning function is calculating wrongly")
                     assert False 
   
##################For testing the generated global nodes and numbering the element with corresponding global number##################
def test_globnodesandnumber():
    print ("Unit_Test_4a")  
    br=3;L=3;l=1
    #Breadth of the plate(Br)(Horizontal length)-->Corresponds to rows-->x-direction
    #Length of the plate(L)(---->Vertical length)-->Corresponds to columns-->y=direction
    nor= int((br/l));         #Number of elements in a meshed plate in X-direction
    noc= int((L/l));          #Number of elements in a meshed plate in Y-direction
    aDOFsiz=0;usize=0;
    globU=[];
    globnodesextra=[];
    new=[]
    globnodes,aDOFsiz,usize,globnodesextra,globU,u,aDOF,cornernodenum,allboundarynodenum,ordinarynodenum=naming_umatrix(noc,nor,l,new,aDOFsiz,usize,globnodesextra,globU)    
    expected_val=([[0,0],[1,0],[2,0],[3,0],[0,1],[1,1],[2,1],[3,1],[0,2],[1,2],[2,2],[3,2],[0,3],[1,3],[2,3],[3,3]])
    if np.array_equal(globnodes,expected_val):
        print("The generated global coordinates for the given model is correct")
        assert True
    else:
        print("The generated global coordinates for the given model is not correct")
        assert False
    print ("\n")
    print("Unit_Test_4b")
    count=0;elementc=[]
    elementswglobnode,globnumber=glob_number(noc,nor,count,elementc,u,aDOF) 
    expected_val1=([[4,0,1,5],[5,1,2,6],[6,2,3,7],[8,4,5,9],[9,5,6,10],[10,6,7,11],[12,8,9,13],[13,9,10,14],[14,10,11,15]])   
    if np.array_equal(globnumber,expected_val1):
        print("The generated global numbers for the respective elements are correct")
        assert True
    else:
        print("The generated global numbers for the respective elements are not correct")
        assert False
