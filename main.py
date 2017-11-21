# -*- coding: utf-8 -*-
def f(x1, x2, x3, x4):
    """Retorna o valor da função f(x) no ponto (x1 x2 x3 x4)"""
    return -30*x1 - 10*x1*x2 - 2*x1*x3 - 3*x1*x4 - 10*x2 - 10*x2*x3 - 10*x2*x4 - 40*x3 - x3*x4 - 12*x4


def df(x1, x2, x3, x4):
    """Retorna a derivada de f(x) no ponto (x1 x2 x3 x4)"""
    return [-30-10*x2-2*x3-3*x4, -10*x1-10-10*x3-10*x4, -2*x1-10*x2-40-x4, -3*x1-10*x2-x3-12]

print df(0,0,0,0)