#Quiz 5

print("Devin Fan")


# PROBLEM 0
import numpy as np
import numpy.linalg as la
import scipy.integrate as si
import matplotlib.pyplot as plt


# PROBLEM 1


def my_erf(x):
    x=np.linspace(0,x,100)
    y=np.exp(-x^2)
    F=2/np.sqrt(np.pi)*si.simps(x,y)
    return F

my_erf(1/4)


# PROBLEM 2

r=np.array([0,.8,1.2,1.4,2.0,3,3.4,3.6,4.0,5,5.5,6.4])
r=r*10
p=np.array([1.8,1.79,1.77,1.7,1.67,1.56,1.49,1.45,1.23,1.18,1.15,1.13])

dv=4*np.pi*r**2
m=si.simps(p,dv)
print(m)


# PROBLEM 3

m=-0.5

def fx(x,v,t): return v
def fv(x,v,t): return -x+(m*(1-x**2)*v)

h=0.05
t=np.arange(0,40+h,h)

x=np.zeros(t.shape)
v=np.zeros(t.shape)

x[0]=2
v[0]=0

for i in range(len(t)-1):
    
    kx1=fx(x[i],v[i],t[i])
    kv1=fv(x[i],v[i],t[i])
    
    kx2=fx(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    kv2=fv(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    
    kx3=fx(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    kv3=fv(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    
    kx4=fx(x[i]+h*kx3,v[i]+h*kv3,t[i])
    kv4=fv(x[i]+h*kx3,v[i]+h*kv3,t[i])
    
    kx=(kx1+2*kx2+2*kx3+kx4)/6
    kv=(kv1+2*kv2+2*kv3+kv4)/6
    
    x[i+1]=x[i]+h*kx
    v[i+1]=v[i]+h*kv

plt.plot(t,x,'b-')
plt.title('damped harmonic oscillator')
plt.ylabel('x')
plt.xlabel('t')
plt.show()

m=5

def fx(x,v,t): return v
def fv(x,v,t): return -x+(m*(1-x**2)*v)

h=0.05
t=np.arange(0,40+h,h)

x=np.zeros(t.shape)
v=np.zeros(t.shape)

x[0]=0
v[0]=2

for i in range(len(t)-1):
    
    kx1=fx(x[i],v[i],t[i])
    kv1=fv(x[i],v[i],t[i])
    
    kx2=fx(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    kv2=fv(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    
    kx3=fx(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    kv3=fv(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    
    kx4=fx(x[i]+h*kx3,v[i]+h*kv3,t[i])
    kv4=fv(x[i]+h*kx3,v[i]+h*kv3,t[i])
    
    kx=(kx1+2*kx2+2*kx3+kx4)/6
    kv=(kv1+2*kv2+2*kv3+kv4)/6
    
    x[i+1]=x[i]+h*kx
    v[i+1]=v[i]+h*kv

w=plt.plot(t,x,'b-')
plt.title('damped harmonic oscillator')
plt.ylabel('x')
plt.xlabel('t')
plt.show()

m=5

def fx(x,v,t): return v
def fv(x,v,t): return -x+(m*(1-x**2)*v)

h=0.05
t=np.arange(0,40+h,h)

x=np.zeros(t.shape)
v=np.zeros(t.shape)

x[0]=0
v[0]=2

for i in range(len(t)-1):
    
    kx1=fx(x[i],v[i],t[i])
    kv1=fv(x[i],v[i],t[i])
    
    kx2=fx(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    kv2=fv(x[i]+h*kx1/2,v[i]+h*kv1/2,t[i]+h/2)
    
    kx3=fx(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    kv3=fv(x[i]+h*kx2/2,v[i]+h*kv2/2,t[i]+h/2)
    
    kx4=fx(x[i]+h*kx3,v[i]+h*kv3,t[i])
    kv4=fv(x[i]+h*kx3,v[i]+h*kv3,t[i])
    
    kx=(kx1+2*kx2+2*kx3+kx4)/6
    kv=(kv1+2*kv2+2*kv3+kv4)/6
    
    x[i+1]=x[i]+h*kx
    v[i+1]=v[i]+h*kv

w=plt.plot(x,t,'b-')
plt.title('damped harmonic oscillator')
plt.ylabel('x')
plt.xlabel('t')
plt.show()


