#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import math
import array

L = 100/0.05292 #atm
N = 300
#V_0 = 0
delta_x = L/N
m = 0.067 #atm
h_cross = 1

E_meV = np.linspace(0,35, 1000)
x = np.zeros(N)
psi = np.zeros(N)
V = np.zeros(N)

for i in range (N):
    x[i] = delta_x * i
    
def psipsi(E,V):
    #psi[0] = 0
    psi[1] = 1
    #psi[N] = 0
    
    for i in range(1, N-1):
        psi[i+1] = -(2*m *(E - V[i]) *delta_x**2 *psi[i]) - psi[i-1] + 2*psi[i]
        
    C = 0
    for i in range(N):
        C += delta_x * psi[i]**2
    
    for i in range(N):
        psi[i] = psi[i]/math.sqrt(C)
        
    return psi[N-1]

psi_values = []
E_atm = E_meV / 27211.6 #atm

for E in (E_atm):
    
    psi_values.append(psipsi(E,V))

    


# In[3]:


print(psi_values)
print(E_atm)


# In[2]:


#E_meV = E_range * 27211.6 # meV
plt.figure
plt.plot(E_meV, psi_values)
size = (len(psi_values))


# In[3]:


#x = np.linspace(0, L,size)

E_zero_points_index = np.argmin(np.abs(np.array(psi_values))) 
E_zero_points = E_atm[E_zero_points_index]

E_up_5 = E_zero_points * 1.05
E_down_5 = E_zero_points * 0.95

psi_e_0 = psipsi(E_zero_points,V)
psi_e_pl_5 = psipsi(E_up_5,V)
psi_e_min_5 = psipsi(E_down_5,V)

print(E_zero_points)
print(E_up_5)
print(E_down_5)


# In[ ]:





# In[4]:


#bisection
th = 0.000001
def bisection(a,b):
    
    #f_middle = 0
    f_a = psipsi(a,V)
    f_b = psipsi(b,V)
    
    if f_a*f_b>0: 
        print("guwno guwno zjedz je ruwno kutasiarzu")
    
    for _ in range (N): 
        middle = (a+b)/2
        f_middle = psipsi(middle,V)
        
        if (b-a)/2 <th:
            return middle
        
        if f_middle*f_a < 0:
            b = middle
            f_b = f_middle
        else:
            a = middle
            f_a = f_middle
    return (a+b)/2
            


# In[5]:


#def find_0s():
#    zero_points = []

zero_points = []
for i in range(len(E_atm) -1 ):
    if psipsi(E_atm[i],V)*psipsi(E_atm[i + 1],V) < 0:
        zero_point = bisection(E_atm[i], E_atm[i+1])
        zero_points.append(zero_point)
        
zero_points_eV = np.array(zero_points) * 27211.6
print(zero_points_eV)


# In[6]:


def zero_teor(n):
    E_teor = np.zeros(n)
    const = math.pi**2/(2*m*L**2)
    for i in range(1,n+1):
        E_teor[i-1] = i*i*const*27211.6
    return E_teor

E_teor = zero_teor(7)
print(E_teor)


# In[7]:


W_val = np.linspace(0,1000,200) #meV
zero_points_W = []
E_atm_2 = np.linspace(0-50,35, 1000)/27211.6
psi = np.zeros(N)

def potential(W): 
    V1 = np.zeros(N)
    V1[N//2]=-W/27211.6
    return V1

for W in (W_val):
    V = potential(W)
    zero_points_2 = []
    
    for i in range(len(E_atm_2) - 1):
        if psipsi(E_atm_2[i],V)*psipsi(E_atm_2[i + 1],V) < 0:
            zero_point_2 = bisection(E_atm_2[i], E_atm_2[i+1])
            #print(zero_point_2)
            zero_points_2.append(zero_point_2)
        if len(zero_points_2)>= 7:
            break
            
    zero_points_eV_2 = np.array(zero_points_2) * 27211.6
    zero_points_W.append(zero_points_eV_2)
    
for i in range(100):  
    print(zero_points_W[i])


# In[8]:


plt.figure()
for i, W in enumerate(W_val):
    for j in range(len(zero_points_W[i])):
        plt.scatter(W, zero_points_W[i][j], color="red")
    
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




