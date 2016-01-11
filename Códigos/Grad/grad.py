#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Bibliotecas importadas
'''
import numpy as np
from scipy import linalg as la

# Dados de entrada
eps = 1e-3  # Tolerância
v0 = v = np.array([5.,1.25])  # Chute inicial
A = [np.array([[3.,2.],[2.,6.]]), np.array([[1.,0.],[0.,1.]])]  # Matrizes A
b = [np.array([-2.,-8.]), np.array([-2.,2.])]  # Vetores b
k = 30  # Número limites de iterações
mat_num = 0  # Index para selecionar qual matriz será usada

# Cálculo dos autovalores e autovetores
evals_A1, evecs_A1 = la.eig(A[0])
evals_A2, evecs_A2 = la.eig(A[1])
print "Autovalores de A1:\n", evals_A1
print "Autovetores de A1:\n", evecs_A1
print "Autovalores de A2:\n", evals_A2
print "Autovetores de A2:\n", evecs_A2

# Método dos gradientes
print "-----Métodos dos gradientes-----"

x = v0[0]  # Vetor para armazenar as coordenadas das iterações em x
y = v0[1]  # Vetor para armazenar as coordenadas das iterações em y
r = np.dot(A[mat_num],v) + b[mat_num]
for i in range(1,k):
  print "\nIteração: %i" % (i)
  v_ant = v
  print "v_ant = ", v_ant
  tmin = np.dot(r,r)/np.dot(np.dot(A[mat_num],r),r)
  print "tmin = ", tmin
  v = v - tmin*r
  x = np.append(x,v[0])
  y = np.append(y,v[1])
  print "v = ", v
  r = r - tmin*np.dot(A[mat_num],r)
  print "r = ", r
  if ((la.norm(r, np.inf)<= eps) or (la.norm((v-v_ant), np.inf)/la.norm(v, np.inf)<= eps)):
    vsol = v
    print "\n-----------------------------------"
    print "Solução: ", vsol
    print "Chute inicial: ", v0
    print "Tolerância utilizada = %.3f\nNúmero de iterações = %d" % (eps,i)
    print "Os valores dos candidatos em cada iteração foram:"
    print "x:\n", x
    print "y:\n", y
    break

# Método dos gradientes conjugados
print "\n-----Métodos dos gradientes conjugados-----"

flag = 0
x = v0[0]  # Vetor para armazenar as coordenadas das iterações em x
y = v0[1]  # Vetor para armazenar as coordenadas das iterações em y
print "\nIteração: 1"
r = np.dot(A[mat_num],v0) + b[mat_num]
r_ant2 = r
p = -1*r
print "p = ", p
q = np.dot(r,r)/np.dot(np.dot(A[mat_num],r),r)
print "q = ", q
v = v0 + q*p
x = np.append(x,v[0])
y = np.append(y,v[1])
print "v = ", v
r = r + q*np.dot(A[mat_num],p)
r_ant1 = r
print "r = ", r
if ((la.norm((v-v0), np.inf)/la.norm(v, np.inf)<= eps) or (la.norm(r,np.inf)<=eps)):
  vsol = v
  flag = 1
  print "\n-----------------------------------"
  print "Solução: ", vsol
  print "Chute inicial: ", v0
  print "Tolerância utilizada = %.3f\nNúmero de iterações = %d" % (eps,1)
  print "Os valores dos candidatos em cada iteração foram:"
  print "x:\n", x
  print "y:\n", y

if flag == 0:
  for j in range(2,k):
    print "\nIteração: %i" % (j)
    v_ant = v
    alpha_ant = np.dot(r_ant1,r_ant1)/np.dot(r_ant2,r_ant2)
    p = -r + alpha_ant*p
    print "p = ", p
    q = np.dot(r,r)/np.dot(np.dot(A[mat_num],p),p)
    print "q = ", q
    v = v + q*p
    x = np.append(x,v[0])
    y = np.append(y,v[1])
    print "v = ", v
    r_ant2 = r_ant1
    r = r + q*np.dot(A[mat_num],p)
    r_ant1 = r
    print "r = ", r
    if ((la.norm((v-v_ant), np.inf)/la.norm(v, np.inf)<= eps)):
      vsol = v
      print "\n-----------------------------------"
      print "Solução: ", vsol
      print "Chute inicial: ", v0
      print "Tolerância utilizada = %.3f\nNúmero de iterações = %d" % (eps,j)
      print "Os valores dos candidatos em cada iteração foram:"
      print "x:\n", x
      print "y:\n", y
      break
