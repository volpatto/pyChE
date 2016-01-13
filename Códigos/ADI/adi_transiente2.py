#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from scipy import linalg as la
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt

def sol(x,y):
    '''
    Solução exata
    '''
    return (np.sin(np.pi*x)*np.sin(np.pi*y))/(2.*np.pi**(2.))

def f(i,j,h):
    '''
    Termo de geração
    '''
    return (np.sin(np.pi*h*i)*np.sin(np.pi*h*j))

itmax = 5000
N = 30
Ne = 200
h = 1./(N)
tf = 0.5
Dt = 0.01
tole = 10**(-10.)
A = np.zeros((N,N))
C = np.zeros((N,N))
u = np.zeros((N,N))
utemp_i = np.zeros((N,N))
utemp_j = np.zeros((N,N))
b = np.zeros(N)
d = np.zeros(N)

# Parâmetros iterativos transiente
sigx = Dt/(2.*h**2.)
sigy = Dt/(2.*h**2.)

x_pe = np.linspace(0, 1., N-1)
y_pe = np.linspace(0, 1., N-1)
Xe,Ye = np.meshgrid(x_pe, y_pe)
sol_exata = sol(Xe,Ye)

# Inicialização
t = np.array([])
t = np.append(t,0.)

# Inicializando o plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# Definindo o grid da aproximação
x_p = np.linspace(0, 1., N)
y_p = np.linspace(0, 1., N)
X,Y = np.meshgrid(x_p, y_p)
Z = u

#Plotando
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)

fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,y)')
ax.set_title('t = %f' % t[0])
ax.set_zlim([-0.01,0.05])
plt.savefig('adi_surf{0:04d}.png'.format(0))
#plt.show()
plt.clf()

plt.pcolor(X,Y,Z,antialiased=True, vmin=0.0, vmax=0.05)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('t = %f' % t[0])
plt.savefig('adi_cont{0:04d}.png'.format(0))
#plt.show()
plt.clf()

it = 0
while t[it] <= tf:
    it += 1
    if (it)*Dt <= tf:
        t = np.append(t,((it)*Dt))
        for j in range(1,N-1):
            for i in range(N-1):
                if i == 0:
                    A[i,i] = 1.
                if ((i<N-1) and i>0):
                    A[i,i-1] = -sigx
                    A[i,i] = 1 + 2.*sigx
                    A[i,i+1] = -sigx
                    b[i] = sigy*u[i,j-1]+(1-2.*sigy)*u[i,j]+sigy*u[i,j+1]+(Dt/2.)*f(i,j,h)
                if i == (N-1):
                    A[i,i] = 1.

            utemp_i[:,j] = la.solve(A,b)

        u = utemp_i

        for i in range(1,N-1):
            for j in range(N-1):
                if j == 0:
                    C[j,j] = 1.
                if ((j<N-1) and j>0):
                    C[j,j-1] = -sigy
                    C[j,j] = 1 + 2.*sigy
                    C[j,j+1] = -sigy
                    d[j] = sigx*u[i-1,j]+(1-2.*sigx)*u[i,j]+sigx*u[i+1,j]+(Dt/2.)*f(i,j,h)
                if j == (N-1):
                    C[j,j] = 1.

            utemp_j[i,:] = la.solve(C,d)

        u = utemp_j

        x_pe = np.linspace(0, 1., N)
        y_pe = np.linspace(0, 1., N)
        Xe,Ye = np.meshgrid(x_pe, y_pe)
        sol_exata = sol(Xe,Ye)

        print u
        apx = la.norm((u-sol_exata),np.inf)/la.norm((sol_exata),np.inf)
        print apx

        # Inicializando o plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        #Plotando
        Z = u
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=True,vmax=0.05)

        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u(x,y)')
        ax.set_title('t = %f' % t[it])
        if (np.abs(u.max()) and np.abs(u.min()))<0.05:
            ax.set_zlim([-0.01,0.05])
        else:
            ax.autoscale()
        plt.savefig('adi_surf{0:04d}.png'.format(it))
        y_pexata = np.linspace(0, 1., Ne)
        plt.clf()

        x_pexata = np.linspace(0, 1., Ne)
        sol_exata = sol(Xe,Ye)
        Xe,Ye = np.meshgrid(x_pexata, y_pexata)

        plt.pcolor(X,Y,Z,antialiased=True, vmin=0.0, vmax=0.05)
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('t = %f' % t[it])
        plt.savefig('adi_cont{0:04d}.png'.format(it))
        plt.clf()

    else:
        break
