#!/usr/bin/env python
# coding: utf-8

# In[42]:


import numpy as np
import matplotlib.pyplot as plt
import math
import array

L = 100/0.05292 #atm
N = 100
W = 0
delta_x = L/N
m = 0.067 #atm
h_cross = 1
alpha = delta_x**2*m/h_cross**2 *0.95

E_meV = np.linspace(0,35, 1000)
V = np.zeros(N+1)
x = np.zeros(N+1)

for i in range (N+1):
    x[i] = delta_x * i
    

def normalization(psi):
    I = 0
    for i in range(N+1):
        I += delta_x * abs(psi[i])**2
    psi /= math.sqrt(I)  
    return psi
#    for i in range(N):
#        psi[i] = psi[i]/math.sqrt(I)    
    
def hamiltonian(psi):
    h_psi = np.zeros_like(psi)
    for i in range(1,N):
        h_psi[i] = -h_cross**2/(2*m) *(psi[i+1]+psi[i-1]-2*psi[i])/delta_x**2 + V[i]*psi[i]
    return h_psi
    
def itm(psi, psi_prim, minus=None):
        
    psi = normalization(psi)
    h_psi = hamiltonian(psi)
    
    E_old = 0
    for i in range(N+1):
        E_old += psi[i]*h_psi[i]*delta_x
    E_old = E_old * 27211.6 #meV
    
    for i in range(1,N):
        psi_prim[i] = psi[i] - alpha*h_psi[i]
#----------------------------
    if minus is not None: 
        c=0
        for i in range(N+1):
            c += minus[i]*psi_prim[i]*delta_x
        psi_prim = psi_prim - c*minus       
#----------------------------        
    psi_prim = normalization(psi_prim)
    h_psi_prim = hamiltonian(psi_prim)
    
    E_new = 0
    for i in range(N+1):
        E_new += psi_prim[i]*h_psi_prim[i]*delta_x
    E_new = E_new * 27211.6
    
    psi,psi_prim = psi_prim, psi
    
        
    delta_E = abs(E_new - E_old) 
    
    
    return psi, delta_E, E_new
        
        
        
        
        
        



# In[48]:


psi_prim = np.zeros(N+1)
psi_prim2 = np.zeros(N+1)
#h_psi = np.zeros(N+1)
psi1 = np.random.uniform(-1,1,N+1)
psi1[0] = 0
psi1[N] = 0
psi2 = np.random.uniform(-1,1,N+1)
psi2[0] = 0
psi2[N] = 0
E = []
k=0
c=0
while 1>0:
    psi1,delta_E, E_new = itm(psi1,psi_prim)
    k+=1
    E.append(E_new)
    if delta_E < 1e-6:
        #print(k, delta_E)
        break
while 1>0:
    psi2,delta_E2,E_new2 = itm(psi2,psi_prim2,psi1)
    if delta_E2 < 1e-6:
        break
print(E)       
plt.plot(x,psi)


# In[40]:


plt.plot(np.linspace(1,len(E),len(E)), E)
plt.xscale("log")
plt.yscale("log")


# In[49]:


plt.plot(x, psi2)


# In[ ]:




