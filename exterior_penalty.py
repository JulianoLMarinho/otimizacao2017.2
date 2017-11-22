# -*- coding: utf-8 -*-
import sympy as S
from scipy import arange

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

for i in arange(-0.5,2,0.5):
    for j in arange(-0.5, 2,0.5):
        for k in arange(-0.5, 2,0.5):
            for l in arange(-0.5, 2, 0.5):
                print "X=("+str(i)+","+str(j)+","+str(k)+","+str(l)+")" + "\tf(X) = " +str(f.subs(w,i).subs(x,j).subs(y,k).subs(z,l)) + "\th(X) = " + str(h.subs(w,i).subs(x,j).subs(y,k).subs(z,l))