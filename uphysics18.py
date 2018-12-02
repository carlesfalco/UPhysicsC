# Python modules
from numpy import linspace
import random as random
from scipy.integrate import odeint
from math import sqrt, pi, log
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.size': 13})
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['tahoma']

# Equations of motion spherical sail
def light_spherical_sail(y,t,Frad,Lc,R,m,I):
    x, vx, thetay, omegay = y
    dydt = [vx, Frad(t)/m*thetay-1/2/m*Frad(t)*(x+Lc*thetay)/R,
            omegay, -1/2*Frad(t)/I*(x+Lc*thetay)*(Lc/R-1)]
    return dydt

# Equations of motion flat sail
def light_flat_sail(y,t,Frad,Lc,m,I):
    x, vx, thetay, omegay = y
    dydt = [vx, Frad(t)/m*thetay,
            omegay, 1/2*Frad(t)/I*(x+Lc*thetay)]
    return dydt
	
	# Constants
a = 200 #cm
R = 1000 #cm
Lc = 1000 #cm
m1 = 1/2 #g
m2 = 1/2 #g
m = (m1+m2) #g
L = 2000 #cm
I = m1*Lc**2+m2*(L-Lc)**2 #g*cm^2
F_rad_0 = 10**7 #cm*g/s^2
tim = 10*60 #s
P = 50*10**9 #W
c = 3*10**8 #m/s
eta = 1 #reflectance
def F_rad(x): #Force from the laser beam
    if x<=10*60:
        return F_rad_0  #+normal(F_rad_0/100,F_rad_0/100) #adding noise
    if x>10*60:
        return 0
		
y0 = [0.001,0,0.001,0] #Initial cond set 1
y1 = [0.0001,0,0.0001,0] #Initial cond set 2
yr = [0.0000001,0,0.00001,0] #Initial cond set 3

sol1 = odeint(light_spherical_sail, y0, t, args = (F_rad,Lc,R,m,I))
sol2 = odeint(light_spherical_sail, y0, t2, args = (F_rad,Lc,R,m,I))
sol12 = odeint(light_spherical_sail, y1, t, args = (F_rad,Lc,R,m,I))
sol22 = odeint(light_spherical_sail, y1, t2, args = (F_rad,Lc,R,m,I))
# solving for different conditions

#Plotting results from ODEs
fig,ax = plt.subplots(5,1,figsize=(7,12))

ax[0].plot(t3,z)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('z (km)')

ax[1].plot(t,sol1[:,0],label=r'$\delta_{x} = 0.001$')
ax[1].legend()
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('x (cm)')
ax[1].set_xlim(0,tim/2)
ax[1].set_ylim(-1.5,3.5)

ax[2].plot(t2,sol2[:,0])
ax[2].set_xlabel('Time (s)')
ax[2].set_ylabel('x (cm)')
ax[2].set_xlim(0,tim/300)
ax[2].set_ylim(-1,3)

ax[3].plot(t,sol12[:,0],label=r'$\delta_{x} = 0.0001$',color='orange')
ax[3].legend()
ax[3].set_xlabel('Time (s)')
ax[3].set_ylabel('x (cm)')
ax[3].set_xlim(0,tim/2)
ax[3].set_ylim(-0.15,0.35)

ax[4].plot(t2,sol22[:,0],color='orange')
ax[4].set_xlabel('Time (s)')
ax[4].set_ylabel('x (cm)')
ax[4].set_xlim(0,tim/300)
ax[4].set_ylim(-0.1,0.3)

plt.savefig('dxbis.eps')
plt.show()

#Solving for flat sail
sol_flat = odeint(light_flat_sail,y1,t4,args = (F_rad,Lc,m,I))
fig,ax = plt.subplots(figsize=(7,4.5))
ax.plot(t4,sol_flat[:,1])
plt.show()