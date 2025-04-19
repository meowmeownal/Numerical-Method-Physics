#!/usr/bin/env python
# coding: utf-8

# In[79]:


import numpy as np
import matplotlib.pyplot as plt
import math
import array
import cmath

delta_x = 0.5 / 0.05292
m = 0.067 #m_0
E = 7 / 27211.6
E_2 = np.linspace(1,50, 20000) / 27211.6
#E_2 = list(E_2)

#q =  math.sqrt(2*E*m)

barrier = 10 / 27211.6

psi = np.zeros(340, dtype=complex)
V = np.zeros(340)
x = np.zeros(340)

for i in range (340):
    x[i] = delta_x * i

for i in range(0,100):
    V[i] = 0
for i in range(100,120):
    V[i] = barrier
for i in range(120,220):
    V[i] = 0
for i in range (220,240):
    V[i] = barrier
for i in range (240,340):
    V[i] = 0
    
def wavef(E):
    #q = math.sqrt(2*E*m)
    psi[-1] = 1
    psi[-2] = cmath.exp(-1j*q*delta_x)
    for i in range (338,0,-1):
        psi[i-1]= -2*m*(E - V[i])*delta_x**2*psi[i] - psi[i+1] + 2*psi[i]
    return psi
    
def stale(psi):
    A = ( psi[0]* cmath.exp(1j*q*x[0]) - psi[1]*cmath.exp(1j*q*x[1]) ) / ( (cmath.exp(1j*q*x[0]))**2 - (cmath.exp(1j*q*x[1]))**2 )
    B = -1*((-psi[1]*cmath.exp(1j*q*x[0]) + cmath.exp(1j*q*x[1])*psi[0]) * cmath.exp(1j*q*x[1]+1j*q*x[0]))/((cmath.exp(1j*q*x[0]))**2-(cmath.exp(1j*q*x[1]))**2)
    R = abs(B)**2/abs(A)**2
    T = 1/abs(A)**2
    return A,B,R,T
    


# In[80]:


q = math.sqrt(2*E*m)
psi= wavef(E)
#print(E_2 *27211.6)


# In[81]:


psi_ab = np.zeros_like(psi)
A,B,R,T = stale(psi)
for i in range (340):
    psi_ab[i] = A*cmath.exp(1j*q*x[i]) + B*cmath.exp(-1j*q*x[i])


# In[82]:


fig, ax1 = plt.subplots()
ax1.plot(x*0.05292,np.abs(psi)**2)
ax2 = ax1.twinx()
ax2.step(x*0.05292, V*27211.6, color = "orange")
ax1.set_xlabel("x nm")
ax1.set_ylabel("nowa fala grubson")
ax2.set_ylabel("potencjał niewyczerpany chyba w dna on byl mi dany")
ax1.plot(x*0.05292,np.abs(psi_ab)**2, color = "pink", linestyle='--')


# In[83]:


psi_TR = []
A2 =[]
B2 = []
R2 = []
T2 = []

for  E in E_2:
    q = math.sqrt(2*E*m)
    psi_TR.append(wavef(E))
    A,B,R,T = stale(psi)
    T2.append(T)
    R2.append(R)


#print(T2)


# In[84]:


maximaT = []
k =[]
T2 = np.array(T2)
index =[]

for i in range (len(T2)-1):
    k.append( (T2[i+1] - T2[i])/(E_2[i+1] - E_2[i]) )
    
for i in range (len(k)-1):
    if k[i]*k[i+1] < 0:
        if T2[i]>0.98:
            index.append(i)
print(index)


# In[85]:


E_2_meV = E_2*27211.6
plt.plot(E_2_meV, T2)
plt.plot(E_2*27211.6, R2)
plt.plot(E_2[index][:4]*27211.6, T2[index][:4], "xr")


# In[ ]:



        


# In[ ]:




