#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
(Portuguese only)
About:
O presente codigo foi desenvolvido para resolver numericamente os modelos de Thomas e de Bohart-Adams de adsorcao em leitos fixos
dinamicos. Utiliza-se o metodo de diferencas finitas com o esquema upwind. Eh necessario observar o cumprimento da condicao CFL do
metodo, tendo em vista a natureza convectiva do problema. O codigo nao sera executado caso tal condicao nao seja suprida. No presente
caso, a condicao eh a seguinte:
						|u_0|*dt/dh < 1
				
Author: Diego Tavares Volpatto
Contact: dtvolpatto@gmail.com
'''

# Bibliotecas importadas

import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt
import sys
import time
from scipy import interpolate as intp

def fread(file_name, var_str, var_kind):
    '''
    Esta função lê o arquivo file_name (string), busca a string var_str (string) e lê o valor associado e o retorna. A variável var_kind
    especifica o tipo da variável a ser lida.

    Exemplo:
    No arquivo, tem-se:

    var_qualquer{
        23.5
    }

    Fazendo:
    var = fread('nomearquivo.txt','var_qualquer'), var recebe 23.5.

    Autor: Diego T. Volpatto
    E-mail: dtvolpatto@gmail.com
    '''

    f = open(file_name, 'r')
    read_data = f.readlines()
    if (var_kind == 'float'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                var = float(read_data[line+1])
                return var
    if (var_kind == 'int'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                var = int(read_data[line+1])
                return var
    if (var_kind == 'bool'):
        for line in range(len(read_data)):
            if ('#' in read_data[line]): 
                continue
            if (var_str in read_data[line]):
                var = bool(read_data[line+1])
                return var

def main():

# Leitura dos parametros do modelo

    arq_entrada = 'input.txt'

    BA_model = fread(arq_entrada,'BA_model', 'int')
    if (BA_model == None): 
        sys.exit("Erro leitura externa BA_model")
    else:
        BA_model = bool(BA_model)

    Thomas_model = fread(arq_entrada,'Thomas_model', 'int')
    if (Thomas_model == None): 
        sys.exit("Erro leitura externa Thomas_model")
    else:
        Thomas_model = bool(Thomas_model)

    C_0 = fread(arq_entrada,'C_0', 'float')
    if (C_0 == None): 
        sys.exit("Erro leitura externa C_0")

    d = fread(arq_entrada,'id', 'float')
    if (d == None): 
        sys.exit("Erro leitura externa d")

    h = fread(arq_entrada,'height', 'float')
    if (h == None): 
        sys.exit("Erro leitura externa h")

    Q = fread(arq_entrada,'Q', 'float')
    if (Q == None): 
        sys.exit("Erro leitura externa Q")
    else:
        u_0 = Q/(np.pi*(d/2.0)**2.0)	# cm/min

    eps = fread(arq_entrada,'eps', 'float')
    if (eps == None): 
        sys.exit("Erro leitura externa eps")

    rho_s = fread(arq_entrada,'rho_s', 'float')
    if (rho_s == None): 
        sys.exit("Erro leitura externa rho_s")

    if (Thomas_model == True):

        b = fread(arq_entrada,'b', 'float')
        if (b == None): 
            sys.exit("Erro leitura externa b")

        q_max = fread(arq_entrada,'q_max', 'float')
        if (q_max == None): 
            sys.exit("Erro leitura externa q_max")

        k_a1 = fread(arq_entrada,'k_a1', 'float')
        if (k_a1 == None): 
            sys.exit("Erro leitura externa k Thomas")

    if (BA_model == True):
        
        k_a2 = fread(arq_entrada,'k_a2', 'float')
        if (k_a2 == None): 
            sys.exit("Erro leitura externa k BA")

        q_e = fread(arq_entrada,'q_e', 'float')
        if (q_e == None): 
            sys.exit("Erro leitura externa q_e")

# Parametros do metodo numerico

    m = fread(arq_entrada,'nump', 'int')
    if (m == None): 
        sys.exit("Erro leitura externa nump")

    dt = fread(arq_entrada,'dt', 'float')
    if (dt == None): 
        sys.exit("Erro leitura externa dt")

    tf = fread(arq_entrada,'t_total', 'float')
    if (tf == None): 
        sys.exit("Erro leitura externa tf")

    Np = m + 1	# Numero de pontos discretos
    dh = h/m	# Passo espacial
    tole = 1.0e-6	# Tolerancia do criterio do ponto de ruptura

    if (Thomas_model == True):
        C_T = np.zeros(Np)	# Vetor de solucao da C
        q_T = np.zeros(Np)	# Vetor de solucao da q
        CnewT = np.copy(C_T)
        qnewT = np.copy(q_T)
    if (BA_model == True):
        C_BA = np.zeros(Np)	# Vetor de solucao da C
        q_BA = np.zeros(Np)	# Vetor de solucao da q
        CnewBA = np.copy(C_BA)
        qnewBA = np.copy(q_BA)

# Dados experimentais

    if (C_0 == 0.6):
        exptime, expdata = np.loadtxt('C06.txt', unpack=True)
    elif (C_0 == 1.2):
        exptime, expdata = np.loadtxt('C12.txt', unpack=True)
    else:
        exptime, expdata = np.loadtxt('C16.txt', unpack=True)

# Inicialização

    t = np.array([])
    t = np.append(t,0.)
    if (Thomas_model == True):
        Ch_T = np.array([])
        Ch_T = np.append(Ch_T,0.)
        ruptureT = False
    if (BA_model == True):
        Ch_BA = np.array([])
        Ch_BA = np.append(Ch_BA,0.)
        ruptureBA = False

# Teste de condicao CFL de estabilidade

    cfl = np.abs(u_0)*(dt/dh)

# Inicializacao de variaveis graficas

    t_print = np.loadtxt('tprint.txt', unpack=True)

    z = np.linspace(0.0,h,Np)

    if (Thomas_model == True):
        fig1 = plt.figure('CC0T')
        ax1 = plt.subplot(111)

        fig2 = plt.figure('qT')
        ax2 = plt.subplot(111)

    if (BA_model == True):
        fig3 = plt.figure('CC0BA')
        ax3 = plt.subplot(111)

        fig4 = plt.figure('qBA')
        ax4 = plt.subplot(111)

# Loop de resolucao

    it = 0
    while t[it] <= tf:
        if (cfl >= 1.0):
                print u"Condição CFL não cumprida. Não há estabilidade."
                sys.exit("Erro: CFL")
                #break
        it += 1
        if (it)*dt <= tf:
            t = np.append(t,((it)*dt))
            if (t[it] != 0.0):
                    if (Thomas_model == True):
                        CnewT[0] = C_0
                    if (BA_model == True):
                        CnewBA[0] = C_0
            else:
                    if (Thomas_model == True):
                        CnewT[0] = C_0
                    if (BA_model == True):
                        CnewBA[0] = C_0
            for i in range(1,m+1):			
                    if (Thomas_model == True):
                        # Modelo de Thomas
                        CnewT[i] = C_T[i] - dt*u_0*(C_T[i]-C_T[i-1])/dh - dt*(rho_s/eps)*k_a1*((q_max-q_T[i])*C_T[i]-(1.0/b)*q_T[i])
                        qnewT[i] = q_T[i] + dt*k_a1*(q_max-q_T[i])*CnewT[i] - dt*(k_a1/b)*q_T[i]
                    if (BA_model == True):
                        # Modelo de Bohart-Adams
                        CnewBA[i] = C_BA[i] - dt*u_0*(C_BA[i]-C_BA[i-1])/dh - dt*(rho_s/eps)*(k_a2*(q_e-q_BA[i])*C_BA[i])
                        qnewBA[i] = q_BA[i] + dt*k_a2*(q_e-q_BA[i])*CnewBA[i]
                    if (i == 1):
                        if (Thomas_model == True):
                            qnewT[i-1] = q_T[i-1] + dt*k_a1*(q_max-q_T[i-1])*CnewT[i-1] - dt*(k_a1/b)*q_T[i-1]
                        if (BA_model == True):
                            qnewBA[i-1] = q_BA[i-1] + dt*k_a2*(q_e-q_BA[i-1])*CnewBA[i-1]
            if (it*dt in t_print):
                if (Thomas_model == True):
                    ax1.plot(z,CnewT,linewidth=2.5,label=("%d min" % (int(it*dt))),figure=fig1)
                    ax2.plot(z,qnewT,linewidth=2.5,label=("%d min" % (int(it*dt))),figure=fig2)
                if (BA_model == True):
                    ax3.plot(z,CnewBA,linewidth=2.5,label=("%d min" % (int(it*dt))),figure=fig3)
                    ax4.plot(z,qnewBA,linewidth=2.5,label=("%d min" % (int(it*dt))),figure=fig4)
                    
            if (Thomas_model == True):
                C_T = np.copy(CnewT)
                q_T = np.copy(qnewT)
                Ch_T = np.append(Ch_T,(C_T[-1]/C_T[0]))
                resCh_T = np.abs(Ch_T[-2] - Ch_T[-1])/np.abs(Ch_T[-1])
                if (resCh_T < tole and ruptureT == False):
                        rpoint_T = t[-1]
                        ruptureT = True
            if (BA_model == True):
                C_BA = np.copy(CnewBA)
                q_BA = np.copy(qnewBA)
                Ch_BA = np.append(Ch_BA,(C_BA[-1]/C_BA[0]))
                resCh_BA = np.abs(Ch_BA[-2] - Ch_BA[-1])/np.abs(Ch_BA[-1])
                if (resCh_BA < tole and ruptureBA == False):
                        rpoint_BA = t[-1]
                        ruptureBA = True
            print u"Passo = %d; Tempo = %f min " % (it, (it*dt))
        else:
            break

    print u"\n--------------------------------\n"
    print u"Condição CFL = %f" % cfl
    print u"Concentração na entrada da coluna: C_0 = %f" % C_0

    if (Thomas_model == True):
        print u"--- Parâmetros do modelo de Thomas ---"
        print u"b = %f" % b
        print u"q_max = %f" % q_max
        print u"k_a1 = %f" % k_a1
        if (ruptureT == True): 
            print u"Tempo de ruptura Thomas = %f" % rpoint_T

    if (BA_model == True):
        print u"--- Parâmetros do modelo de Bohart-Adams ---"
        print u"q_e = %f" % q_e
        print u"k_a2 = %f" % k_a2
        if (ruptureBA == True): 
            print u"Tempo de ruptura Bohart-Adams = %f" % rpoint_BA

# Confeccao dos graficos

    if (Thomas_model == True):

# Plot do perfil de concentracao no leito (Thomas)

        ax1.set_xlabel(r'$z\,(cm)$',fontsize=16)
        ax1.set_ylabel(r'$C\,\left(\frac{mmol}{L}\right)$',fontsize=18)
        box1 = ax1.get_position()
        ax1.set_position([0.1*box1.x0+box1.x0, 0.1*box1.y0 + box1.y0, box1.width * 0.75, box1.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.grid(True)

# Plot do perfil de concentração adsorvida no leito (Thomas)

        ax2.set_xlabel(r'$z\,(cm)$',fontsize=16)
        ax2.set_ylabel(r'$\overline{q}\,\left(\frac{mmol}{g}\right)$',fontsize=18)
        box2 = ax2.get_position()
        ax2.set_position([0.1*box2.x0+box2.x0, 0.1*box2.y0 + box2.y0, box2.width * 0.75, box2.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax2.grid(True)

    if (BA_model == True):

# Plot do perfil de concentracao no leito (Bohart-Adams)

        ax3.set_xlabel(r'$z\,(cm)$',fontsize=16)
        ax3.set_ylabel(r'$C\,\left(\frac{mmol}{L}\right)$',fontsize=18)
        box3 = ax3.get_position()
        ax3.set_position([0.1*box3.x0+box3.x0, 0.1*box3.y0 + box3.y0, box3.width * 0.75, box3.height])
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax3.grid(True)

# Plot do perfil de concentração adsorvida no leito (Bohart-Adams)

        ax4.set_xlabel(r'$z\,(cm)$',fontsize=16)
        ax4.set_ylabel(r'$\overline{q}\,\left(\frac{mmol}{g}\right)$',fontsize=18)
        box4 = ax4.get_position()
        ax4.set_position([0.1*box4.x0+box4.x0, 0.1*box4.y0 + box4.y0, box4.width * 0.75, box4.height])
        ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax4.grid(True)

# Plot das curvas de rupturas dos modelos

    fig = plt.figure("rupture")
    ax = plt.subplot(111)

    if (Thomas_model == True):
        apSpl_T = intp.splrep(t,Ch_T,s=0)
        ap_T = intp.splev(exptime,apSpl_T,der=0)
        ax.plot(t,Ch_T,'-.',label="Thomas",linewidth=2.5,color='blue')
        ax.plot(exptime,ap_T,'x',label="Thomas (previsto)",linewidth=2.5)
        np.savetxt(("prev_T_%d.dat" % int(C_0*10)),ap_T,fmt='%.15e')
        np.savetxt(("Ch_T_%d.dat" % int(C_0*10)),Ch_T,fmt='%.15e')
    if ((C_0 == 0.6) or (C_0 == 1.2) or (C_0 == 1.6)):
        ax.plot(exptime,expdata,'o',label="Valores experimentais",linewidth=2.5,color='blue')
    if (BA_model == True):
        apSpl_BA = intp.splrep(t,Ch_BA,s=0)
        ap_BA = intp.splev(exptime,apSpl_BA,der=0)
        ax.plot(t,Ch_BA,'-.',label="Bohart-Adams",linewidth=2.5,color='red')
        ax.plot(exptime,ap_BA,'x',label="Bohart-Adams (previsto)",linewidth=2.5)
        np.savetxt(("prev_BA_%d.dat" % int(C_0*10)),ap_BA,fmt='%.15e')
        np.savetxt(("Ch_BA_%d.dat" % int(C_0*10)),Ch_BA,fmt='%.15e')
    np.savetxt(("t_%d.dat" % int(C_0*10)),t,fmt='%.15e')
    ax.set_xlabel(r'$t\,(min)$',fontsize=16)
    ax.set_ylabel(r'$C/C_0$',fontsize=18)
    box = ax.get_position()
    ax.set_position([0.1*box.x0+box.x0, 0.1*box.y0 + box.y0, box.width * 0.75, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#ax.set_ylim([-0.1,Ch_BA.max()+0.1*Ch_BA.max()])
    ax.set_ylim([-0.1,1.0])
    plt.title('Curva de ruptura')
    plt.grid(True)
    plt.show()

main()
