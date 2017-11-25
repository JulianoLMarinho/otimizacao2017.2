# -*- coding: utf-8 -*-
import operator
import sympy as S
from scipy import arange
import numpy as np

w, x, y, z = S.symbols('w,x,y,z')
f = -30 * w - 10 * w * x - 2 * w * y - 3 * w * z - 10 * x - 10 * x * y - 10 * x * z - 40 * y - y * z - 12 * z
# 33w + 14x + 47y + 11z <= 59
# Xi >= 0, min Xi, 0
# Xi <= 1 Max Xi-1, 0

# Criar uma nova funcao f adicionando penalidade

p = 1 # fator a ser alterado posteriormente
g1 = 33*w + 14*x + 47*y + 11*z - 59
g2 = w
g3 = w-1
g4 = x
g5 = x-1
g6 = y
g7 = y-1
g8 = z
g9 = z-1

h = f + p* (S.Max(0,g1)**2 + S.Min(0,g2)**2+S.Max(0,g3)**2+S.Min(0,g4)**2+S.Max(0,g5)**2+S.Min(0,g6)**2+S.Max(0,g7)**2+S.Min(0,g8)**2+S.Max(0,g9)**2)

print "Função inicial: f(X) = " + str(f) + "\nFunção com penalidade exterior: h(X) = " + str(h)

def multEscalarVetores(v1,v2):
    return np.dot(np.array(v1),np.array(v2))

def metodoGradiente(func, x0):
    k = 0
    xk = x0
    t = S.symbols('t')
    gradiente = [S.diff(func,w).subs(w, xk[0]).subs(x,xk[1]).subs(y, xk[2]).subs(z,xk[3]),S.diff(func,x).subs(w, xk[0]).subs(x,xk[1]).subs(y, xk[2]).subs(z,xk[3]),S.diff(func,y).subs(w, xk[0]).subs(x,xk[1]).subs(y, xk[2]).subs(z,xk[3]),S.diff(func,z).subs(w, xk[0]).subs(x,xk[1]).subs(y, xk[2]).subs(z,xk[3])]
    print gradiente
    while(all(v==0 for v in gradiente)==False and k<=0):
        td = [m*(-1*t) for m in gradiente]
        print "vetor td"
        print td
        xktd = [xk[i]+td[i] for i in range(len(xk))]
        print "vetor xk+tkdk"
        print xktd
        tmp = h.subs(w, xktd[0]).subs(x, xktd[1]).subs(y, xktd[2]).subs(z, xktd[3])
        print "equacao penalizada substituida"
        print tmp.subs(t, 3)
        # temos a funcao ja em fi(t).
        k=k+1


    # p1, p2, p3, p4 = ponto
    # # calcula o gradiente e pega seu valor no ponto fornecido
    # grad1=S.diff(func,w).subs(w, p1).subs(x, p2).subs(y, p3).subs(z, p4)
    # grad2=S.diff(func,x).subs(w, p1).subs(x, p2).subs(y, p3).subs(z, p4)
    # grad3=S.diff(func,y).subs(w, p1).subs(x, p2).subs(y, p3).subs(z, p4)
    # grad4=S.diff(func,z).subs(w, p1).subs(x, p2).subs(y, p3).subs(z, p4)
    # d = dinicial
    # # gera uma lista com os valores do gradiente
    # point = [grad1,grad2,grad3,grad4]
    # #mult escalar de vetores
    # temp = multEscalarVetores(point, d)
    # print temp

#buscaIntervalo(h, [1,0,1,0], [0.5,0.5,0.5,0.5])

#metodoGradiente(h, [2,2,2,2])

for i in arange(-0.5,2,0.5):
    for j in arange(-0.5, 2,0.5):
        for k in arange(-0.5, 2,0.5):
            for l in arange(-0.5, 2, 0.5):
                print "X=("+str(i)+","+str(j)+","+str(k)+","+str(l)+")" + "\tf(X) = " +str(f.subs(w,i).subs(x,j).subs(y,k).subs(z,l)) + "\th(X) = " + str(h.subs(w,i).subs(x,j).subs(y,k).subs(z,l))
