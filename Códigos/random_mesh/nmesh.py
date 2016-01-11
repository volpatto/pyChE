#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy import linalg as la
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt
from numpy import random as r

def f(x,y):
    return (np.sin(np.pi*x)*np.sin(np.pi*y))

def sol(x,y):
    return (np.sin(np.pi*x)*np.sin(np.pi*y))/(2.*np.pi**(2.))

N = 30
h = 1./N
itmax = 5000
tole = 10.**(-10.)
w = 2. - 2.*np.pi*h
u = np.zeros((N,N))
x1 = np.linspace(0,1,N)
y1 = np.linspace(0,1,N)
XX,YY = np.meshgrid(x1,y1)
print XX
print YY
plt.scatter(XX,YY)
plt.savefig('nmesh_reg{0:02d}.png'.format(N))
plt.show()

for i in range(1,N-1):
    for j in range(1,N-1):
        if r.randint(0,2) == 1:
            s = 1.
        else:
            s = -1.
        phi = s*r.random()
        XX[i,j] = XX[i,j] + phi*(h/3.)

print XX

for j in range(1,N-1):
    for i in range(1,N-1):
        if r.randint(0,2) == 1:
            s = 1.
        else:
            s = -1.
        phi = s*r.random()
        YY[i,j] = YY[i,j] + phi*(h/3.)

print YY

plt.scatter(XX,YY)
plt.savefig('nmesh_nonu{0:02d}.png'.format(N))
plt.show()

# Solver SOR
for it in range(itmax+1):
    print u"iteração = %d" % it
    u_ant = np.copy(u)
    for i in range(1,N-1):
        for j in range(1,N-1):
            coord0 = np.array([XX[i,j-1],YY[i,j-1]])
            coord1 = np.array([XX[i,j],YY[i,j]])
            coord2 = np.array([XX[i,j+1],YY[i,j+1]])
            coord3 = np.array([XX[i-1,j],YY[i-1,j]])
            coord4 = np.array([XX[i+1,j],YY[i+1,j]])
            h1 = la.norm(coord1-coord2)
            h2 = la.norm(coord1-coord0)
            h3 = la.norm(coord1-coord3)
            h4 = la.norm(coord1-coord4)
            fator = (1./(h1*h2)+1./(h3*h4))**(-1.)
            u[i,j] = fator*((1./(h1+h2))*(u[i,j-1]/h2+u[i,j+1]/h1)+\
            (1./(h3+h4))*(u[i-1,j]/h4+u[i+1,j]/h3)+(1./2.)*f(coord1[0],coord1[1]))
            u[i,j] = w*u[i,j]+(1.-w)*u_ant[i,j]
    conv = la.norm((u-u_ant),np.inf)/la.norm(u,np.inf)
    if conv <= tole:
        print "conv = ", conv
        break

print u
sol_exata = sol(XX,YY)
apx = la.norm((u-sol_exata),np.inf)/la.norm((sol_exata),np.inf)
print u"Erro na norma do máx = %f" % apx
print u"núm de iterações = %d" % it

# Inicializando o plot
fig = plt.figure()
ax = fig.gca(projection='3d')

Z = u

#Plotando
surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
plt.savefig('nmesh_surf{0:02d}.png'.format(N))
plt.show()

plt.pcolor(XX,YY,Z)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('nmesh_cont{0:02d}.png'.format(N))
plt.show()

x2 = np.linspace(0.,1.,200)
y2 = np.linspace(0.,1.,200)
Xe,Ye = np.meshgrid(x2,y2)
Ze = sol(Xe,Ye)
plt.pcolor(Xe,Ye,Ze)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('nmesh_exact.png')
plt.show()
