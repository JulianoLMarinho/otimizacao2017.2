# -*- coding: utf-8 -*-
import sympy as S
# f(x1, x2, x3, x4) = -30x1 - 10x1x2 - 2x1x3 - 3x1x4 - 10x2 - 10x2x3 -10x2x4 - 40x3 - x3x4 - 12x4
# x1=w
# x2=x
# x3=y
# x4=z

w, x, y, z = S.symbols('w,x,y,z')


f = -30 * w - 10 * w * x - 2 * w * y - 3 * w * z - 10 * x - 10 * x * y - 10 * x * z - 40 * y - y * z - 12 * z

gw = S.diff(f, w)
gx = S.diff(f, x)
gy = S.diff(f, y)
gz = S.diff(f, z)

criticos = S.solve([gw,gx,gy,gz],[w,x,y,z])
mh = S.hessian(f, [w,x,y,z])

# autovalores agora
autovalores=S.Matrix.berkowitz_minors(mh)

print "Seja a função f"
S.pprint(f)
print "\nCalculando o gradiente, temos: "
print [gw,gx,gy,gz]
print "\nEncontrando os pontos criticos desse gradiente: "
print criticos
print "\nCalculamos então a matriz hessiana:"
S.pprint(mh)
print "\nPor fim, encontramos os autovalores dessa matriz"
print autovalores
