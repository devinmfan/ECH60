#!/usr/bin/env python
# coding: utf-8

# In[80]:


import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

def rsquared(x,y,yfit):
    ymean=np.mean(y)
    sstot=sum((y-ymean)**2)
    ssres=sum((y-yfit)**2)
    rs=1-ssres/sstot
    return rs

def bisection(f,a,b,tol,max_iter):
    for i in range(max_iter):
        x_est=(a+b)/2    
        if f(a)*f(x_est)<0:
            b=x_est
        else:
            a=x_est
        if np.abs(f(x_est))<tol:
            break
    if i==max_iter-1: 
        print('No solution after %d iterations' %max_iter)
        return
    else:
        return x_est
        print('The numerical solution is', x_est, 'after', i, 'iterations')


# In[81]:


k1=0.05
k2=0.15
k3=0.03
w1=0.5
w2=0.3
w3=0.2
hin=1.0
hout=0.8
tin=200
tout=20

Linear_equation_matrix=[[((k1/w1)+hin),-k1/w1,0,0],[-(k1/w1),((k2/w2)+(k1/w1)),-k2/w2,0],[0,-(k2/w2),((k3/w3)+(k2/w2)),-k3/w3],[0,0,-(k3/w3),(hout+(k3/w3))]]
print(Linear_equation_matrix)

Constant_matrix=[[hin*tin],[0],[0],[hout*tout]]
print(Constant_matrix)

T0,T1,T2,T3=la.solve(Linear_equation_matrix,Constant_matrix)
print(T0)
print(T1)
print(T2)
print(T3)


# In[83]:


k1_2=0.05
k2_2=0.15
k3_2=0.03
w2_2=0.3
w3_2=0.2
hin_2=1.0
hout_2=0.8
tin_2=200
tout_2=20

for w1_3 in range(1,21,1):
  w1_2=w1_3/10
  Linear_equation_matrix2=[[((k1_2/w1_2)+hin_2),-(k1_2/w1_2),0,0],[-(k1_2/w1_2),((k2_2/w2_2)+(k1_2/w1_2)),-(k2_2/w2_2),0],[0,-(k2_2/w2_2),((k3_2/w3_2)+(k2_2/w2_2)),-(k3_2/w3_2)],[0,0,-(k3_2/w3_2),(hout_2+(k3_2/w3_2))]]
  Constant_matrix2=[[hin_2*tin_2],[0],[0],[hout_2*tout_2]]
  T0,T1,T2,T3=la.solve(Linear_equation_matrix2,Constant_matrix2)
  print(T3)


# In[86]:


ydata=np.array([37.41935484,35.08379888,33.30049261,31.89427313,30.75697211,29.81818182,29.03010033,28.35913313,27.78097983,27.27762803,26.83544304,26.44391408,26.09480813,25.78158458,25.49898167,25.24271845,25.00927644,24.79573712,24.59965928,24.41898527])
print(ydata)


# In[87]:


x=np.arange(0.1,2.1,0.1)
y=ydata[0:]

plt.figure(figsize=(6,6))

plt.plot(x,y,'ro')
plt.xlabel('Thickness of Layer #1 (cm)')
plt.ylabel('Temperature of outside edge, T3 (Degrees Celcius)')

m,n=np.polyfit(np.log(x), y, 1)
plt.plot(x,m*np.log(x)+n)

print('coefficient of determination is %.4f' %rsquared(x,y,m*np.log(x)+n))


# In[75]:


f=lambda p: m*np.log(p)+(n-30)
fmin=lambda p: m*np.log(p)+(n-29.5)
fmax=lambda p: m*np.log(p)+(n-30.5)
root=bisection(f,0.000001,1,0.000001,100)
max_width=bisection(fmin,0.000001,1,0.000001,100)
min_width=bisection(fmax,0.000001,1,0.000001,100)
max= '{0:.4g}'.format(max_width)
min= '{0:.4g}'.format(min_width)


# In[77]:


print('For the instrument to work, the minimum width necessary is', min, 'cm and the maximum is', max, 'cm.')

