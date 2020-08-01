# Prandtl's-lifting-line-solution
Discretized the wing into N segments. Prandtl's lifting-line theory was used to model the system.  Included more segments for each harmonic, so as to resolve the issue of illconditioned matrix. Least square method was used to solve the matrix equation.  Bound circulation was evaluated for both elliptical and tapered wings. 
"""
Created on Wed Jul  1 07:50:36 2020

@author: Souridas A
"""
"""
================================================
Given data:
    Zero lift angle of attack = 0 degree
    lift curve slope of airfoil = 2*pi/radian
    Free stream velocity = 10m/s
    mid-span chord length = 0.1m
    Wing span = 1m
    Angle of attack = 5 degree
================================================
"""
"""
In the first part of the code,provide N=40,M=40 as mentioned in the problem.This creates unexpected results.The model created is underfitted.
In the second part, use N=40,M=100. To solve the problem of underfitting, here the number of observations(positions) are increased.
"""
import matplotlib.pyplot as plt
import numpy as np
import math
#N denotes the number of frequencies
N=input("Enter number of frequencies required : ")
N=int(N)
#M denotes the number of sections
M=input("Enter number of sections required : ")
M=int(M)
pi=math.pi
alpha=5*pi/180  #Angle of attack in radians
#Calculating the midpoints  and corresponding theta value of each section
positions=np.linspace(-.5,.5,M+1)
midpoints=np.empty([M,1])
thetas=np.empty([M,1])
for i in range(M):
    midpoints[i]=(positions[i]+positions[i+1])/2
    thetas[i]=math.acos(-2*midpoints[i]/1)
#Elliptical Chord distribution fuction
def chord_ellipse_distr(y):
    return 0.1*math.sqrt(1-4*y**2)
#Function which creates the bound circulation function.It's parameter is the chord distribution function
def bound_circ(distr_func):
    M_matrix=np.empty([M,N])
    for i in range(M):
        for j in range(N):
            M_matrix[i][j]=(2*1+.5*2*pi*distr_func(midpoints[i])*(j+1)/math.sin(thetas[i]))*math.sin(thetas[i]*(j+1))
    #RHS Matrix creation
    RHS=np.empty([M,1])
    for i in range(M):
        RHS[i]=.5*2*pi*distr_func(midpoints[i])*alpha
    # Creating A matrix(least square method is employed)
    A=np.dot(np.dot(np.linalg.inv(np.dot(M_matrix.transpose(),M_matrix)),M_matrix.transpose()),RHS)
    # finding bound circulation
    Bound_circ=np.empty([M,1])
    for j in range(M):
           sum=0
           for i in range(N):
               sum=sum+2*1*10*A[i]*math.sin(thetas[j]*(i+1))
           Bound_circ[j]=sum
    return Bound_circ,A


# Comparing numerical and analytical solution for bound circulation by plotting(ELLIPTICAL CHORD DISTRIBUTION)
plt.title("NUMERICAL AND ANALYTIC SOLUTION FOR BOUND CIRCULATION (ELLIPTICAL CHORD)")
plt.xlabel("Span Positions (in meters)")
plt.ylabel("Bound Circulation")
plt.scatter(midpoints,bound_circ(chord_ellipse_distr)[0],c='Red')
gamma_0=pi*alpha/(1+.1*pi/2)
def analytic_sol(y):
    return gamma_0*math.sqrt(1-4*y**2)
analytic_points=np.empty([M,1])
for i in range(M):
    analytic_points[i]=analytic_sol(midpoints[i])
plt.plot(midpoints,analytic_points)
plt.legend(['ANALYTIC SOLUTION','NUMERICAL SOLUTION'])
plt.show()
#Bar graph for first 10 A values
A_10=['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10']   
A_list=[]
for i in range(10):
    A_list.append(bound_circ(chord_ellipse_distr)[1][i][0])
plt.ylabel("First 10 fourier coefficents")
plt.bar(A_10,A_list,color='red')
plt.show()


#SOLUTION TO QUESTION 2
#Creating chord distribution functions for each tapered wing
def tapered_1(y):
        return .1
def tapered_2(y):
        return .1*(1-abs(y))
def tapered_3(y):
        return .1-.2*abs(2*y)/3
             
#Plotting numerical solution for bound circulation (Tapered wing)
plt.title("NUMERICAL SOLUTION FOR BOUND CIRCULATION (TAPERED CHORD)")
plt.xlabel("Span Positions (in meters)")
plt.ylabel("Bound Circulation")
plt.plot(midpoints,bound_circ(tapered_1)[0],c='Red')
plt.plot(midpoints,bound_circ(tapered_2)[0],c='Blue')
plt.plot(midpoints,bound_circ(tapered_3)[0],c='Orange')
plt.legend(['taper ratio=1:1','taper ratio=2:1','taper ratio=3:1'])
plt.show()
#plotting first 10 A values for 3 wings
A1_list=[]    
A2_list=[]
A3_list=[]

for i in range(10):
    A1_list.append(bound_circ(tapered_1)[1][i][0])   
    A2_list.append(bound_circ(tapered_2)[1][i][0]) 
    A3_list.append(bound_circ(tapered_3)[1][i][0]) 
    
ind=np.arange(10)
width=.35 
plt.bar(ind,A1_list,width,label='taper ratio=1:1')
plt.bar(ind+width,A2_list,width,label='taper ratio=2:1')
plt.bar(ind+2*width,A3_list,width,label='taper ratio=3:1')
plt.ylabel('First 10 fourier coefficients')
plt.xticks(ind+width,('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10'))
plt.legend(loc='best')
plt.show()
