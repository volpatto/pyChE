#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Bibliotecas importadas
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

def g1(v1,v2):
  '''
  Essa função representa, computacionalmente, o funcional
  da equação matricial do item 1
  '''
  return ((3./2.)*v1**2. + 3.*v2**2. + 2.*v1*v2 - 2.*v1 - 8.*v2)

def g2(v1,v2):
  '''
  Essa função representa, computacionalmente, o funcional
  da equação matricial do item 2
  '''
  return (v1**2.)/2. - 2.*v1 + (v2**2.)/2. + 2.*v2

# Criando os pontos e os grids de plot
x1 = np.linspace(-4.0, 3.0, 200)
y1 = np.linspace(-2.0, 5.0, 200)
x2 = np.linspace(-2.0, 6.0, 200)
y2 = np.linspace(-6.0, 2.0, 200)
X1, Y1 = np.meshgrid(x1, y1)
X2, Y2 = np.meshgrid(x2, y2)

# Configurando e plotando as curvas de nível de q1
CS1 = plt.contour(X1, Y1, g1(X1,Y1), 15, colors='k')
'''
# Pontos do Gradiente
x = np.array([ 1.,          1.15846995, -0.06052963,
              -0.03277463, -0.24627463, -0.24141352,
              -0.27880668, -0.27795529, -0.28450446,
              -0.28435535, -0.28550239])
y = np.array([-1.,          0.90163934,  1.00322264,
              1.33628264,  1.3540743,   1.41240764,
              1.41552373,  1.42574044,  1.42628621,
              1.4280756,   1.42817119])
plt.plot(x,y,color='k')
'''
'''
# Pontos do Gradiente Conjugado
x = np.array([ 1.,          1.15846995, -0.28571429,
               -0.28571429])
y = np.array([-1.,          0.90163934,  1.42857143,
                1.42857143])
plt.plot(x,y,color='k')
'''
plt.xlabel('x')
plt.ylabel('y')
plt.clabel(CS1, inline=1, fontsize=10)
plt.title(u'Curvas de nível de q(v) para 1')
xevec1 = np.array([-2.,1.])
base = np.array([-0.28571429,1.42857143])
xevec2 = np.array([1.,2.])
q1 = quiver(base[0],base[1],xevec1[0],xevec1[1],scale=8,width=0.005)
plt.text(0.25,2.0,'7')
q2 = quiver(base[0],base[1],xevec2[0],xevec2[1],scale=8,width=0.005)
plt.text(-1.2,1.5,'2')
plt.show()

# Configurando e plotando as curvas de nível de q2
CS2 = plt.contour(X2, Y2, g2(X2,Y2), 15, colors='k')
'''
# Pontos do Gradiente e do GC
x = np.array([ 5.,  2.])
y = np.array([1.25, -2.])
plt.plot(x,y,color='k')
'''
plt.xlabel('x')
plt.ylabel('y')
plt.clabel(CS2, inline=1, fontsize=10)
plt.title(u'Curvas de nível de q(v) para 2')
xevec1 = np.array([1.,0.])
base = np.array([2.,-2.])
xevec2 = np.array([0.,1.])
q1 = quiver(base[0],base[1],xevec1[0],xevec1[1],scale=8,width=0.005)
plt.text(2.3,-2.3,'1')
q2 = quiver(base[0],base[1],xevec2[0],xevec2[1],scale=8,width=0.005)
plt.text(1.8,-1.6,'1')
plt.show()

# Configurando e plotando o gráfico 3D de g1
fig = plt.figure()
ax = fig.gca(projection='3d')
X3 = np.arange(-4, 3, 0.025)
Y3 = np.arange(-2, 5, 0.025)
X3, Y3 = np.meshgrid(X3, Y3)
Z3 = g1(X3,Y3)
ax.plot_surface(X3, Y3, Z3, rstride=10, cstride=10, color='w')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('q(v) para 1')
plt.show()

# Configurando e plotando o gráfico 3D de g2
fig = plt.figure()
ax = fig.gca(projection='3d')
X4 = np.arange(-2, 6, 0.025)
Y4 = np.arange(-6, 2, 0.025)
X4, Y4 = np.meshgrid(X4, Y4)
Z4 = g2(X4,Y4)
ax.plot_surface(X4, Y4, Z4, rstride=10, cstride=10, color='w')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('q(v) para 2')
plt.show()
