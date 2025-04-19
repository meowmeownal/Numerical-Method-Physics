#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math
import array
import cmath
import numba as nb

m = 0.067 #m_0
h_bar = 1
omega = 5 /27211.6 #meV to atm
dt = 1 #atm
delta_x = 1 / 0.05292 #nm to atm
T = 2*math.pi/omega
x_0 = 30 / 0.05292 #nm to atm

x = np.linspace(-100,100,201) / 0.05292 #nm to atm

@nb.njit
def normalization(psi):
    C = 0
    for i in range(len(psi)):
        C += delta_x * abs(psi[i])**2  
    return C

@nb.njit
def hamiltonian(psi, V):
    h_psi = np.zeros_like(psi)

    for i in range(1,len(psi)-1):
        h_psi[i] = -h_bar**2/(2*m) *(psi[i+1]+psi[i-1]-2*psi[i]) /delta_x**2 + V[i]*psi[i]
    return h_psi
   
def psi_12():
    psi1 = np.zeros(201,dtype=complex)
    psi2 = np.zeros(201,dtype=complex)
    
    psi1[0] = 0
    psi1[200] = 0
    
    for i in range (1,len(psi1)-1):
        psi1[i] = cmath.exp(-m*omega*(x[i]-x_0)**2/2)
        
    #psi1 = psi1/normalization(psi1)
    psi1 /= math.sqrt(normalization(psi1)) 
    psi2 = psi1*cmath.exp(-1j*omega*dt/2)
    
    return psi1, psi2



# In[3]:


V1 = np.zeros(len(x))
for i in range (len(x)):
    V1[i] = m*omega**2*x[i]**2/2

psi_prob = []
psi1,psi2 = psi_12()

for i in range (1,math.floor(10*T)):
    psi1 = psi1 + hamiltonian(psi2,V1)*2*dt/1j/h_bar
    psi1, psi2 = psi2, psi1
    
    if i%500 == 0:
        psi_prob.append(np.abs(psi2)**2)

X,Y = np.meshgrid(x*0.05292, np.linspace(0,len(psi_prob)-1, len(psi_prob))*500)
plt.pcolormesh(X,Y,psi_prob)


# In[29]:


#zad 4
x_expected = []
x_theo =[]

for j in range(len(psi_prob)):
    k = 0
    for i in range (201):
        k+=delta_x*x[i]*psi_prob[j][i]
    x_expected.append(k*0.05292)

size = len(psi_prob)
x_axis = np.linspace(0,size-1,size)

for i in range(len(psi_prob)):
    x_theo.append(x_0*0.05292*math.cos(omega*i*500))

plt.plot(x_axis,x_expected, label = "expected")
plt.plot(x_axis,x_theo, label = "theo")
plt.legend()


# In[3]:


x_0 = 0
psi_prob = []
psi1,psi2 = psi_12()

for i in range (1,math.floor(10*T)):
    psi1 = psi1 + hamiltonian(psi2, V1)*2*dt/1j/h_bar
    psi1, psi2 = psi2, psi1
    
    if i%500 == 0:
        psi_prob.append(np.abs(psi2)**2)

#print(psi_prob)        
plt.imshow(psi_prob)


# In[4]:


x_0= 0
V2 = np.zeros(len(x))
psi_prob = []
psi1,psi2 = psi_12()

for i in range (1,math.floor(10*T)):
    psi1 = psi1 + hamiltonian(psi2, V2)*2*dt/1j/h_bar
    psi1, psi2 = psi2, psi1
    
    if i%500 == 0:
        psi_prob.append(np.abs(psi2)**2)
       
plt.imshow(psi_prob)


# In[ ]:





# In[ ]:




