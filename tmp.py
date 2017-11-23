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

# Sejam os dados comuns do problema
Xk=()
p = 1
err = pow(10,-8)
max_iteracoes=100


def gradiente(funcao):
	"""@return Gradiente da função"""
	return [S.diff(funcao, w),S.diff(funcao,x),S.diff(funcao,y),S.diff(funcao,z)]

def pontos_criticos(funcao):
	"""@return Pontos críticos da Função"""
	return S.solve(gradiente(funcao),[w,x,y,z])

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


def secao_aurea(funcao, eps, ro, d):
	"""@return Ponto t a ser utilizado por métodos de otimização irrestrita"""
	theta1 = (3-sqrt(5))/2
	theta2 = 1-theta1

	# procurar intervalo [a,b]
	a,s,b=0,ro,2*ro
	while(phi(funcao,b,d) < phi(funcao,s,d)):
		a,s,b=s,b,2*b

	#Procuramos um t*
	u = a + theta_1*(b-a)
	v = a + theta_2*(b-a)
	while((b-a)>epslon):
		if(phi(funcao,u,d) < phi(funcao,v,d)):
			b,v = v, u
			u = a + theta_1*(b-a)
		else:
			a,u= u,v
			v = a + theta_2*(b-a)
	t = (u+v)/2
	return t


def metodo_gradiente(funcao, ponto_inicial):
	"""Resolve o problema de minimização da função pelo método do gradiente"""

	# PS: FALTA IMPLEMENTAR AINDA
	print "WARNING: THIS METHOD IS NOT IMPLEMENTED YET!"
	global Xk
	Xk = ponto_inicial
	Xant = (0,0)

	f_i = penalidade_exterior(funcao)
	print "Funcao irrestrita: " + str(f_i)

	while fabs(Xk[0]-Xant[0]) > err and fabs(Xk[1]-Xant[1]) > err and fabs(Xk[2]-Xant[2]) > err and fabs(Xk[3]-Xant[3]) > err and k < max_iteracoes:
		Xant = Xk
		#Calculo o gradiente da funcao irrestrita e aplico o ponto atual
		grad = gradiente(f_i)


#metodo_gradiente(f, (1,1,1,1))