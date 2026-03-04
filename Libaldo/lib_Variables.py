


from sympy import *
 
a, b, c, d, f, g, h, i, j, k, l, m, n, p, q, r, s, t, u, v, w, x, y, z =symbols('a, b, c, d, f, g, h, i, j, k, l, m, n, p, q, r, s, t, u, v, w, x, y, z')  
a0, b0, c0, d0, f0, g0, h0, i0, j0, k0, l0, m0, n0, p0, q0, r0, s0, t0, u0, v0, w0, x0, y0, z0 =symbols('a0, b0, c0, d0, f0, g0, h0, i0, j0, k0, l0, m0, n0, p0, q0, r0, s0, t0, u0, v0, w0, x0, y0, z0')
a1, b1, c1, d1, f1, g1, h1, i1, j1, k1, l1, m1, n1, p1, q1, r1, s1, t1, u1, v1, w1, x1, y1, z1 =symbols('a1, b1, c1, d1, f1, g1, h1, i1, j1, k1, l1, m1, n1, p1, q1, r1, s1, t1, u1, v1, w1, x1, y1, z1')  
a2, b2, c2, d2, f2, g2, h2, i2, j2, k2, l2, m2, n2, p2, q2, r2, s2, t2, u2, v2, w2, x2, y2, z2 =symbols('a2, b2, c2, d2, f2, g2, h2, i2, j2, k2, l2, m2, n2, p2, q2, r2, s2, t2, u2, v2, w2, x2, y2, z2')
a3, b3, c3, d3, f3, g3, h3, i3, j3, k3, l3, m3, n3, p3, q3, r3, s3, t3, u3, v3, w3, x3, y3, z3 =symbols('a3, b3, c3, d3, f3, g3, h3, i3, j3, k3, l3, m3, n3, p3, q3, r3, s3, t3, u3, v3, w3, x3, y3, z3') 

A,B,C,D,L,M,P,R,T,V,W =symbols('A,B,C,D,L,M,P,R,T,V,W')
A1, B1, C1, D1, L1, M1, P1, R1, T1, V1, W1 =symbols('A1, B1, C1, D1, L1, M1, P1, R1, T1, V1, W1')  
A2, B2, C2, D2, L2, M2, P2, R2, T2, V2, W2 =symbols('A2, B2, C2, D2, L2, M2, P2, R2, T2, V2, W2')
A3, B3, C3, D3, L3, M3, P3, R3, T3, V3, W3 =symbols('A3, B3, C3, D3, L3, M3, P3, R3, T3, V3, W3') 

alpha0, alpha1, alpha2, alpha3 = symbols('alpha0, alpha1, alpha2, alpha3')
beta0, beta1, beta2, beta3 = symbols('beta0, beta1, beta2, beta3')
theta0, theta1, theta2, theta3 = symbols('theta0, theta1, theta2, theta3')

dx, dy, dz, dw, dv, du, dt, df, dh, dV, dA =symbols('dx  dy  dz  dw  dv  du  dt  df  dh  dV  dA' ) 
Lvar=[x, y, z, w, v, u, t, f, h, V, A]
dLvar=[dx, dy, dz, dw, dv, du, dt, df, dh, dV, dA]
dalpha,dbeta,dtheta = symbols('d_alpha,d_beta d_theta')

from sympy.abc import * 

def get_diff_var(var):
    # find therespetive differential varible for x,y,z....
    return dLvar[Lvar.index(var)]