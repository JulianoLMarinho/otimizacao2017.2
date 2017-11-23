# -*- coding: utf-8 -*-
def f(x):
    """Retorna o valor da função f(x) no ponto (x1 x2 x3 x4)"""
    return -30*x[0] - 10*x[0]*x[1] - 2*x[0]*x[2] - 3*x[0]*x[3] - 10*x[1] - 10*x[1]*x[2] - 10*x[1]*x[3] - 40*x[2] - x[2]*x[3] - 12*x[3]


def df(x):
    """Retorna a derivada de f(x) no ponto (x1 x2 x3 x4)"""
    return [-30-10*(x[1])-2*(x[2])-3*(x[3]), -10*(x[0])-10-10*x[2]-10*x[3], -2*x[0]-10*x[1]-40-x[3], -3*x[0]-10*x[1]-x[2]-12]


def MEV(X, e):
    Y = []
    for i in range(len(X)):
        Y.append(X[i]*e)
    return Y

def somaVector(X, Y):
    Z = []
    for i in range(len(X)):
        Z.append(X[i]+Y[i])
    return Z

#ponto crítico: x1 = -3.0, x2=-39.0/20, x3=33.0/2 e x4=-29.0/2

def MGradiente(Xo):
    k = 0
    Xk = Xo
    i = True
    t = 1
    while df(Xk)!=[0.0, 0.0, 0.0,0.0]:
        dk = multiEcalarVector(df(Xk),-1)
        m = f(somaVector(Xk, multiEcalarVector(dk, t)))
        n = f(Xk);
        while m>n:
            m = f(somaVector(Xk, multiEcalarVector(dk, t)))
            n = f(Xk);
            t+=0.1
        Xk = somaVector(Xk, multiEcalarVector(dk, t))
        k+=1



x1 = [0.3, 0.1, 0.4, 0.12]
X = [-3.2, -39.0/20, 33.0/2, -29.0/2]
Xo=[1,1,1,1]
Y = [30,10,40,12]

# MGradiente([-3.1, -39.0/20, 33.0/2, -29.0/2])
# print multiEcalarVector(X, 1)
print df(X)
# print f(X)
