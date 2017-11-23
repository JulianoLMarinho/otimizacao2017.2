# -*- coding: utf-8 -*-
import sympy as S
import numpy as np
from scipy import arange

def f(x):
    """Retorna o valor da função f(x) no ponto (x1 x2 x3 x4)"""
    return -30*x[0] - 10*x[0]*x[1] - 2*x[0]*x[2] - 3*x[0]*x[3] - 10*x[1] - 10*x[1]*x[2] - 10*x[1]*x[3] - 40*x[2] - x[2]*x[3] - 12*x[3]


def df(x):
    """Retorna a derivada de f(x) no ponto (x1 x2 x3 x4)"""
    return [-30-10*(x[1])-2*(x[2])-3*(x[3]), -10*(x[0])-10-10*x[2]-10*x[3], -2*x[0]-10*x[1]-40-x[3], -3*x[0]-10*x[1]-x[2]-12]


w, x, y, z = S.symbols('w,x,y,z')
d1, d2, d3, d4, t = S.symbols('d1, d2, d3, d4, t')
f = -30 * w - 10 * w * x - 2 * w * y - 3 * w * z - 10 * x - 10 * x * y - 10 * x * z - 40 * y - y * z - 12 * z

# Criar uma nova funcao f adicionando penalidade

p = 500 # fator a ser alterado posteriormente
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

# print "Função inicial: f(X) = " + str(f) + "\nFunção com penalidade exterior: h(X) = " + str(h)

# for i in arange(-0.5,2,0.5):
#     for j in arange(-0.5, 2,0.5):
#         for k in arange(-0.5, 2,0.5):
#             for l in arange(-0.5, 2, 0.5):
#                 print "X=("+str(i)+","+str(j)+","+str(k)+","+str(l)+")" + "\tf(X) = " +str(f.subs(w,i).subs(x,j).subs(y,k).subs(z,l)) + "\th(X) = " + str(h.subs(w,i).subs(x,j).subs(y,k).subs(z,l))
# print S.diff(h, w)
# print S.diff(h, x)
# print S.diff(h, y)
# print S.diff(h, z)

gf1 = 33*w+t*d1 + 14*x + 47*y + 11*z - 59
gf2 = w+t*d1
gf3 = w+t*d1-1
gf4 = x + t*d2
gf5 = x + t*d2-1
gf6 = y + t*d3
gf7 = y + t*d3-1
gf8 = z + t*d4
gf9 = z + t*d4-1
fi = -30 * (w+t*d1) - 10 * (w+t*d1) * (x + t*d2) - 2 * (w+t*d1) * (y + t*d3) - 3 * (w+t*d1) * (z + t*d4) - 10 * (x + t*d2) - 10 * (x + t*d2) * (y + t*d3) - 10 * (x + t*d2) * (z + t*d4) - 40 * (y + t*d3) - (y + t*d3) * (z + t*d4) - 12 * (z + t*d4)
fip = fi + p* (S.Max(0,g1)**2 + S.Min(0,g2)**2+S.Max(0,g3)**2+S.Min(0,g4)**2+S.Max(0,g5)**2+S.Min(0,g6)**2+S.Max(0,g7)**2+S.Min(0,g8)**2+S.Max(0,g9)**2)
#
#
# k = 0.5*(x-2)**2 + (y-1)**2
#
# kfi = 0.5*((x+t*d1)-2)**2 + ((y+t*d2)-1)**2

# print kfi.subs(x,1).subs(y, 0).subs(d1, 3).subs(d2, 1).subs(t, 0)

# print S.diff(h,w).subs(w,0.5).subs(y,0.5).subs(x,0.5).subs(z,0.5)

# print fip.subs(w,0.5).subs(y,0.5).subs(x,0.5).subs(z,0.5).subs(d1, 1).subs(d2, 1).subs(d3,1).subs(d4, 1)


def tdearmijo(X, d, f, l, N):
    t = 1
    df1 = S.diff(f, w).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3])
    df2 = S.diff(f, x).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3])
    df3 = S.diff(f, y).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3])
    df4 = S.diff(f, z).subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3])
    ddf = df1*d[0] + df2*d[1] + df3*d[2] + df4*d[3]
    f1 = f.subs(w, X[0]+t*d[0]).subs(x, X[1]+t*d[1]).subs(y, X[2]+t*d[2]).subs(z, X[3]+t*d[3])
    f2 = f.subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]) + N*t*ddf
    while f1 > f2:
        t = t*l
        f1 = f.subs(w, X[0] + t * d[0]).subs(x, X[1] + t * d[1]).subs(y, X[2] + t * d[2]).subs(z, X[3] + t * d[3])
        f2 = f.subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]) + N * t * ddf

    return t

def MEV(X, e):
    """multiplica cada elemento de X por e"""
    Y = []
    for i in range(len(X)):
        Y.append(X[i]*e)
    return Y

def SV(X, Y):
    """Soma os vetores X e Y elemento a elemento"""
    Z = []
    for i in range(len(X)):
        Z.append(X[i]+Y[i])
    return Z

def aureo(f):
    O1 = (3.0-5.0**0.5)/2.0
    O2 = 1.0 - O1
    e = 0.0001
    p = 0.5
    a =0.0
    s = p
    b = 2.0*p
    while f.subs(t,b)<f.subs(t,s):
        ft1 = f.subs(t,b)
        ft2 = f.subs(t,s)
        a = s
        s = b
        b = 2.0*b
    u = a + O1*(b-a)
    v = a + O2*(b-a)
    while(b-a)>e:
        if f.subs(t, u)<f.subs(t, v):
            b=v
            v=u
            u=a+O1*(b-a)
        else:
            a = u
            u=v
            v=a+O2*(b-a)
    return (u+v)/2.0


#[0.452753279917195, 0.792590351034160, 0.501284454494597, 0.854621110865347]
def MGradienteA(X, f, armijo):
    k = 0
    Xk = [1.0,1.0,1.0,1.0]
    i = True
    t = 1
    F = f
    dF = [S.diff(f, w), S.diff(f, x), S.diff(f, y), S.diff(f, z)]
    df1 = S.diff(f, w).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
    df2 = S.diff(f, x).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
    df3 = S.diff(f, y).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
    df4 = S.diff(f, z).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
    dFx = [df1, df2, df3, df4]
    e = 0.0001
    while abs(X[0]-Xk[0]) > e:
        dk = MEV(dFx, -1)
        # print "Xk = ", Xk
        if armijo:
            t = tdearmijo(Xk, dk, F, 0.8, 0.25)
            print t
        else:
            H = fip.subs(w, X[0]).subs(x, X[1]).subs(y, X[2]).subs(z, X[3]).subs(d1, dk[0]).subs(d2, dk[1]).subs(d3, dk[2]).subs(d4, dk[3])
            t = aureo(H)
        X = Xk
        Xk = SV(Xk,MEV(dk, t))
        # print "Xk+1 = ", Xk
        # dF = [S.diff(f, w), S.diff(f, x), S.diff(f, y), S.diff(f, z)]
        df1 = S.diff(f, w).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
        df2 = S.diff(f, x).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
        df3 = S.diff(f, y).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
        df4 = S.diff(f, z).subs(w, Xk[0]).subs(x, Xk[1]).subs(y, Xk[2]).subs(z, Xk[3])
        dFx = [df1, df2, df3, df4]
    print Xk
    return Xk



# print tdearmijo([0,0,0,0], [1,1,1,1], h, 0.8, 0.25)
X = [0.0,0.0,0.0,0.0]
Y = [1,1,1,1]
B = [-3.0, -39.0/20, 33.0/2, -29.0/2]

# print MEV(Y, -1)
Xo = MGradienteA(X, h, 1)

# print f.subs(w, Xo[0]).subs(x, Xo[1]).subs(y, Xo[2]).subs(z, Xo[3])

j = (11*t**2)/2 - 5.0*t + 3.0/2

# print j
# print j.subs()
# print aureo(j)