#! /usr/bin/env python
# -*- coding: utf-8 -*-

import methodes
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# ========================================
import sys
reload(sys)
sys.setdefaultencoding('utf8')
# Paramètres du système linéaire, y'' + by' + a^2y = 0, qui est un oscillateur
# harmonique (cas b=0) amorti (b>0) ou amplifié (b<0).
a = 0.9  # fréquence
b = 4. # amortissement
A = np.array([[0.,1.],[-a**2,-b]]) # Matrice du système obtenu par réécriture à
                                   # l'ordre 1 de l'équation
N = 200
T = 20
h = 1/10.
m = [0]*N
n = [0]*N
ha =[0]*N
error = [0]*N
error_r = [0]*N
p1 = [0]*N
p2 = [0]*N
ordr_err = [0]*N
ordr_err_r = [0]*N
C = 10 # constante de complexité
ha[0] = h

# Bornes du plan de phase que l'on va tracer
xmin,xmax, ymin,ymax = -2.,2., -2.,2.

# Pas en x et y utilisé pour tracer les graphes du plan de phase.
hx = (xmax-xmin)/10.
hy = (ymax-ymin)/10.

# Nullclines

# A_11*x + A_12*y = 0
if not A[0,0]==0.:
    y = np.arange(ymin,ymax+hy,hy)
    x = -A[0,1]*y/A[0,0]
else:
    x = np.arange(xmin,xmax+hx,hx)
    y = -A[0,0]*x/A[0,1]

#plt.plot(x,y, linewidth=2)

# A_21*x + A_22*y
if not A[1,1]==0.:
    x = np.arange(xmin,xmax+hx,hx)
    y = -A[1,0]*x/A[1,1]
else:
    y = np.arange(ymin,ymax+hy,hy)
    x = -A[1,1]*y/A[1,0]

#plt.plot(x,y, linewidth=2)

# Champ de vecteur f(x,y):= A*[x,y]^T
y,x = np.mgrid[ymin:ymax+hy:hy, xmin:xmax+hy:hx]
u = A[0,0]*x + A[0,1]*y # Champ des vitesses suivant x
v = A[1,0]*x + A[1,1]*y # Champ des vitesses suivant y
speed = np.sqrt(u*u + v*v) # Champ de vitesses

#plt.streamplot(x,y,u,v, color=speed, linewidth=1)
#plt.axis([xmin,xmax,ymin,ymax])

#plt.colorbar()
#plt.xlabel("y1")
#plt.ylabel("y2)")
#plt.title("Plan de phase système linéaire")
plt.grid(True)

# ========================================
# Partie 2, question 3: code à insérer ici

def Sol(Z, t):
    return np.array([Z[1],-b*Z[1] - a**2 * Z[0]])
def Sol1(Z, t):
    return [Z[1],-b*Z[1] - a**2 * Z[0]]

def erreur_absolue(sol,y0,h):
    x = methodes.euler_explicite(0.,h,N,y0,Sol)
    return np.max( np.abs(x[1]-sol) )
#fonction erreur relative
def erreur_relative(sol,y0,h):
    x = methodes.euler_explicite(0.,h,N,y0,Sol)
    return np.max( np.abs(x[1]-sol) / np.abs(sol) )

t = np.linspace(0,20,200)

for alpha in range(1,4):
 for beta in range(1,7):

  #cas amortie
  if np.power(b,2) > 4*np.power(a,2):
    r_1 = -(1/2) * (b + np.sqrt(np.power(b,2) - 4*np.power(a,2)))
    r_2 = (1/2) * (np.sqrt(np.power(b,2) - 4*np.power(a,2)) - b)
    y0 =  [alpha + beta , alpha*r_1 + beta*r_2]

  #cas harmonique
  elif np.power(b,2) == 4*np.power(a,2):
    y0  =  [alpha + beta , - (b/2) * (alpha + beta)]
    
  #cas amplifié
  else:
    y0 =  [alpha , - (b/2) * alpha + beta * (np.sqrt( 4*np.power(a,2) - np.power(b,2))/2)]
  
  sol_exacte = odeint(Sol1, y0, t )
  sol_appr = methodes.euler_explicite(0,h,N,np.array(y0),Sol)
  err = erreur_absolue(sol_exacte,y0,h)
  err_r = erreur_relative(sol_exacte,y0,h)
  #plt.plot(sol_exacte[:,0], sol_exacte[:,1], 'r')
  plt.plot(t,sol_exacte[:,0])
  #print(err,np.abs(err/sol_exacte))


for k in range(1,20):
    ha[k] = ha[k-1]/10.
    error[k] = erreur_absolue(sol_exacte,y0,ha[k])
    error_r[k] = error[k]/sol_exacte[k]
    p1[k] = (np.log(error[k]) - np.log(error[k-1]))/ float(np.log(ha[k]) - np.log(ha[k-1]))
    p2[k] = np.divide(np.log(error_r[k]) - np.log(error_r[k-1]), np.log(ha[k]) - np.log(ha[k-1]))
    ordr_err[k] = C*np.power(ha[k],p1[k])
    ordr_err_r[k] = C*np.power(ha[k],p2[k])
    print(ha[k] ,error[k],  ordr_err[k],error_r[k] ,  ordr_err_r[k])

#plt.savefig("Err.png")
#plt.show()






# ========================================
