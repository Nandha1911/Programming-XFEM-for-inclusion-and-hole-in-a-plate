import math 

def LevelSet(temp,a,b,r):

  temp1=[];
  for j in range(len(temp)):
        phi=math.sqrt(((temp[j][0]-a)**2)+((temp[j][1]-b)**2))-r; 
        ##Enrichment Function rule:-
        ##If phi value is greater than 0,then the node lies outside the circle and the enrichment 
        ##function will be V(x)=1
        ##If phi value is less than 0,then the node lies inside/on the circle and the enrichment 
        ##function will be V(x)=0
        if phi > 0: #or #phi ==0:
            phi=1;
        elif phi ==0:
            phi = 0;
        else:
            phi=-1;
        temp1.append(phi);                #This 'a' list variable will temporarily rotated inside j loop
                                      #in order to calculate whether the element will be cut by an
                                      #curve or not.....(It proceeds to the if..else..statement)
                                      
  return temp1                               