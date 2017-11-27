# -*- coding: utf-8 -*-
"""
COS360 - Otimização
Trabalho para obtenção de grau na cadeira de Otimização do curso de Engenharia de Computação e Informação - UFRJ - Brasil

@author Lucas Santos de Paula <lucasdepaula@poli.ufrj.br> | Juliano de Lima Marinho <julianomarinho@poli.ufrj.br>

November/2017

"""

from math import *
import sympy as S
from scipy import arange
import numpy as np


# Sejam os símbolos da função w,x,y,z e f a função especificada.
w, x, y, z = S.symbols('w,x,y,z')
f = -30 * w - 10 * w * x - 2 * w * y - 3 * w * z - 10 * x - 10 * x * y - 10 * x * z - 40 * y - y * z - 12 * z
regiao = 33*w+14*x+47*y+11*z-59

# Sejam os dados comuns do problema
Xk=()
p = 5000
err = pow(10,-5)
max_iteracoes=150



class Report():
    """docstring for ClassName"""
    def __init__(self, x0, iteracoes, busca, quant_busca, xotimo, fotimo, erro):
        self.inicial = x0
        self.iteracoes = iteracoes
        self.busca = busca
        self.busca_iter = quant_busca
        self.xotimo = xotimo
        self.fotimo = fotimo
        self.erro = erro
    def __str__(self):
        output = str(self.inicial) + ";" + str(self.iteracoes) + ";" + self.busca + ";" + str(self.busca_iter) + ";" + str(self.xotimo[0]) + ";"+ str(self.xotimo[1]) + ";"+ str(self.xotimo[2]) + ";"+ str(self.xotimo[3]) + ";" + str(self.fotimo) +";" +str(self.erro)
        return str(output)

def moduloV(X):
    return ((X[0]**2)+(X[1]**2)+(X[2]**2)+(X[3]**2))**(1.0/2)

def gradienteV(funcao):
    """@return Gradiente da função"""
    return [S.diff(funcao, w),S.diff(funcao,x),S.diff(funcao,y),S.diff(funcao,z)]

def gradienteE(funcao, X):
    """@return Gradiente da função calculado em X"""
    ret = []
    ret.append(S.diff(funcao, w).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]))
    ret.append(S.diff(funcao, x).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]))
    ret.append(S.diff(funcao, y).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]))
    ret.append(S.diff(funcao, z).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]))
    return ret
#<type 'list'>: [0.590669436471097, 0.560446290980732, 0.612430101224161, 0.545939181145356]
def soma_vetores(X, Y):
    """@return soma os vetores X e Y elemento a elemento"""
    Z = []
    for i in range(len(X)):
        Z.append(X[i]+Y[i])
    return Z

def vetorXescalar(X, esc):
    """@return multiplicação do vetor X pelo escalar esc"""
    temp = []
    for i in range(len(X)):
        temp.append(X[i]*esc)
    return temp

def vetorXvetor(X, Y):
    """@return multiplicação do vetor X pelo vetor Y"""
    temp = 0
    for i in range(len(X)):
        temp+=X[i]*Y[i]
    return temp

def resolve_funcao(funcao, X):
    """@return valor da funcao calculada no ponto X"""
    return funcao.subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3])

def pontos_criticos(funcao):
    """@return Pontos críticos da Função"""
    return S.solve(gradienteV(funcao),[w,x,y,z])

def classificar_pontos_criticos(funcao):
    """"@return True: Imprime se é ponto de sela, mínimo ou máximo. | False: se ocorreu algum erro """
    mh = S.hessian(f, [w,x,y,z])
    autovalores=S.Matrix.berkowitz_minors(mh)
    if all(v>0 for v in autovalores):
        print "Ponto de mínimo"
        return True
    elif all(v<0 for v in autovalores):
        print "Ponto de máximo"
        return True
    else:
        print "Ponto de sela"
        return True
    raise("Não foi possivel classificar os pontos")


def penalidade_exterior(funcao, ro=1):
    """@return Função Penalizada pelo método penalidade externa, transformando um problema restrito em irrestrito."""
    #Sejam as seguintes restricoes
    g1 = 33*w + 14*x + 47*y + 11*z - 59
    g2 = w
    g3 = w-1
    g4 = x
    g5 = x-1
    g6 = y
    g7 = y-1
    g8 = z
    g9 = z-1
    return funcao + ro* (S.Max(0,g1)**2 + S.Min(0,g2)**2+S.Max(0,g3)**2+S.Min(0,g4)**2+S.Max(0,g5)**2+S.Min(0,g6)**2+S.Max(0,g7)**2+S.Min(0,g8)**2+S.Max(0,g9)**2)


def phi(funcao, t, d):
    """@return Valor de phi para o problema de otimização."""
    x1=Xk[0] + t * d[0]
    x2=Xk[1] + t * d[1]
    x3=Xk[2] + t * d[2]
    x4=Xk[3] + t * d[3]
    return funcao.subs(w,x1).subs(x,x2).subs(y,x3).subs(z,x4)


def armijo(funcao, d, gama, eta):
    """@return Ponto t a ser utilizado por métodos de otimização irrestrita"""
    t = 1
    phi1 = phi(funcao, t, d)
    iter=0
    phi2 = np.multiply(vetorXvetor(gradienteE(f, Xk), d), eta*t) + resolve_funcao(funcao, Xk)
    while phi1>phi2:
        iter+=1
        t = t*gama
        phi1 = phi(funcao, t, d)
        phi2 = np.multiply(vetorXvetor(gradienteE(f, Xk), d), eta * t) + resolve_funcao(funcao, Xk)
    return t,iter

def secao_aurea(funcao, erro, ro, d):
    """@return Ponto t a ser utilizado por métodos de otimização irrestrita"""
    theta1 = (3-sqrt(5))/2
    theta2 = 1-theta1

    # procurar intervalo [a,b]
    a,s,b=0,ro,2*ro
    iter = 0
    kkk=0
    while(phi(funcao,b,d) < phi(funcao,s,d)):
        kkk+=1
        a,s,b=s,b,2*b

    #Procuramos um t*
    u = a + theta1*(b-a)
    v = a + theta2*(b-a)
    while((b-a)>erro):
        iter+=1
        if(phi(funcao,u,d) < phi(funcao,v,d)):
            b,v = v, u
            u = a + theta1*(b-a)
        else:
            a,u= u,v
            v = a + theta2*(b-a)
            
    t = (u+v)/2
    return t,iter

def metodo_gradiente(funcao, ponto_inicial, busca):
    """Resolve o problema de minimização da função pelo método do gradiente"""
    # PS: FALTA IMPLEMENTAR AINDA
    global Xk
    totalBusca=0
    Xk = ponto_inicial
    Xant = (50,50, 50, 50)
    k=0
    t = 1
    erro = (0,0,0,0)
    f_i = penalidade_exterior(funcao, p)
    while (fabs(Xk[0]-Xant[0]) > err or fabs(Xk[1]-Xant[1]) > err or fabs(Xk[2]-Xant[2]) > err or fabs(Xk[3]-Xant[3]) > err) and k < max_iteracoes and Xk!=Xant:
        # O critérios de parada utilizados foram X>erro, k máximo
        erro = fabs(moduloV(Xk)-moduloV(Xant))
        Xant = Xk
        print erro
        #Calculo o gradiente da funcao irrestrita e aplico o ponto atual
        d = [i * -1 for i in gradienteV(f_i)]
        d = [i.subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3]) for i in d]

        if busca=='armijo':
            t,q = armijo(f_i, d, 0.8, 0.25)
        elif busca=='aurea':
            t,q = secao_aurea(f_i, 10**(-5), 0.000000001, d)

        totalBusca+=q
        Xk = (Xk[0] + d[0] * t, Xk[1] + d[1] * t, Xk[2] + d[2] * t, Xk[3] + d[3] * t)
        k+=1

    return Report(ponto_inicial,k,busca,totalBusca,Xk,resolve_funcao(f_i, Xk), erro)


arqAu = open("/home/juliano/Projetos/Otimização/result_perto_otimo_aureo.csv","w")#Localização deve ser alterada
arqAr = open("/home/juliano/Projetos/Otimização/result_perto_otimo_armijo.csv","w")
arqAu.write("Xo;Interações.;Busca;Interações da Busca;Opt. X1;Opt. X2;Opt. X3;Opt. X4;Opt. Value; Error\n")
arqAr.write("Xo;Interações.;Busca;Interações da Busca;Opt. X1; Opt. X2;Opt. X3;Opt. X4;Opt. Value; Error\n")
#
Xo=[[0.99,0.98,0.02,0.92],[0.92,0.9,0.0,0.95],[0.9,0.9,0.1,0.9],[1.0,1.0,0.0,1.0], [0.999,0.899,0.02,0.92], [0.93,0.88,0.002,0.89], [0.99,0.99,0.11,0.99], [1.0,0.98,0.1,0.9], [0.965,0.9,0.2,1.0], [0.9,0.93,0.0,0.97], [0.928,0.999,0.0,0.87]]
#
for i in Xo:
    a = str(metodo_gradiente(f, i, 'armijo'))
    b = str(metodo_gradiente(f, i, 'aurea'))
    arqAr.write(a+"\n")
    arqAu.write( b+"\n")
    print a
    print b
# #
arqAr.close()
arqAu.close()

# for i in range(0,11,1):
#     for j in range(0, 11, 1):
#         for k in range(0, 11, 1):
#             for l in range(0, 11, 1):
#                 if(regiao.subs(w,i/10.0).subs(x,j/10.0).subs(y,k/10.0).subs(z,l/10.0))<=0:
#                     arq.write(str([i,j,k,l])+";"
#                               +str(f.subs(w,i/10.0).subs(x,j/10.0).subs(y,k/10.0).subs(z,l/10.0))+";"
#                               +str(regiao.subs(w,i/10.0).subs(x,j/10.0).subs(y,k/10.0).subs(z,l/10.0))+"\n")
#                     print str([i,j,k,l])
# arq.close()
# a = str(metodo_gradiente(f, [10/10.0,10/10.0,10/10.0,10/10.0], 'armijo'))
# b = str(metodo_gradiente(f, [10/10.0,10/10.0,10/10.0,10/10.0], 'aurea'))
# # arq.write(a+"\n")
# # arq.write( b+"\n")
# print a
# print b
# print penalidade_exterior(f,9999999999999999999999999).subs(z,0.296).subs(x,0.292).subs(y,0.704)

# for i in arange(0.00,1.00, 0.01):
#     if (regiao.subs(w, 0.0).subs(x, 0.0).subs(y, i).subs(z, 0.0) <= 0):
#         arq.write(str([1.0, 1.0, i, 1.0]) + ";"
#                   +str(f.subs(w, 0.0).subs(x, 0.0).subs(y, i).subs(z, 0.0))+";"
#                   +str(regiao.subs(w, 0.0).subs(x, 0.0).subs(y, i).subs(z, 0.0))+"\n")
#
# arq.close()
# a = str(metodo_gradiente(f, [1.0,1.0,0.2,1.0], 'armijo'))
# b = str(metodo_gradiente(f, [1.0,1.0,0.2,1.0], 'aurea'))
# # arq.write(a+"\n")
# # arq.write( b+"\n")
# print a
# print b