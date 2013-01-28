#!/usr/bin/python
import scipy, sympy
import itertools

Symbol = sympy.Symbol

I = sympy.I

NS = 4
NC = 3 

eps = [dict(idx = (0,1,2), sign = +1),
       dict(idx = (2,0,1), sign = +1),
       dict(idx = (1,2,0), sign = +1),
       dict(idx = (2,1,0), sign = -1),
       dict(idx = (0,2,1), sign = -1),
       dict(idx = (1,0,2), sign = -1)]

pr = scipy.zeros((NS,NS,NC,NC), dtype=type(Symbol("x")))
pi = scipy.zeros((NS,NS,NC,NC), dtype=type(Symbol("x")))
qr = scipy.zeros((NS,NS,NC,NC), dtype=type(Symbol("x")))
qi = scipy.zeros((NS,NS,NC,NC), dtype=type(Symbol("x")))
for c0 in range(NC):
    for c1 in range(NC):
        for s0 in range(NS):
            for s1 in range(NS):
                pr[s0, s1, c0, c1] = Symbol('A[%d][%d][%d][%d].re' % (s0, s1, c0, c1), real=True)
                pi[s0, s1, c0, c1] = Symbol('A[%d][%d][%d][%d].im' % (s0, s1, c0, c1), real=True)
                                                         
                qr[s0, s1, c0, c1] = Symbol('B[%d][%d][%d][%d].re' % (s0, s1, c0, c1), real=True)
                qi[s0, s1, c0, c1] = Symbol('B[%d][%d][%d][%d].im' % (s0, s1, c0, c1), real=True)

buf = \
    "#ifndef _QCD_PROP_CONTRACT_H\n" + \
    "#define _QCD_PROP_CONTRACT_H 1\n\n"
p = pr + I*pi
q = qr + I*qi
for x0 in range(NS):
    for x1 in range(x0+1, NS):
        x = scipy.zeros([NS,NS,NC,NC], dtype=type(Symbol("x")))
        buf += "__inline__ void\n"
        buf += ("prop_contract_%d%d(qcd_complex_16 C[NS][NS][NC][NC], qcd_complex_16 A[NS][NS][NC][NC], qcd_complex_16 B[NS][NS][NC][NC])\n{\n" 
                 % (x0, x1))
        for col0 in eps:
            a0, b0, c0 = col0["idx"]
            sign0 = col0["sign"] 
            for col1 in eps:
                a1, b1, c1 = col1["idx"]
                sign1 = col1["sign"] 
                for mu in range(NS):
                    for nu in range(NS):
                        for ku in range(NS):
                            idx = [-1]*4
                            idx[x0] = ku
                            idx[x1] = ku
                            idx[idx.index(-1)] = mu
                            idx[idx.index(-1)] = nu
                            csp0 = a1+mu*NC
                            csp1 = a0+nu*NC
                            x[mu, nu, a0, a1] += p[idx[0], idx[1], b0, b1] * q[idx[2], idx[3], c0, c1] * sign0 * sign1

        lines = [x[i,j,k,l].expand().as_real_imag() for i,j,k,l in itertools.product(range(NS), range(NS), range(NC), range(NC))]
        lines = [(str(x[0]), str(x[1])) for x in lines]
        lines = [(" "+x[0] if x[0][0] != "-" else x[0],
                  " "+x[1] if x[1][0] != "-" else x[1]) for x in lines]
        lines = [(str(x[0]).replace(" + ", " \n\t+").replace(" - ", " \n\t-"),
                  str(x[1]).replace(" + ", " \n\t+").replace(" - ", " \n\t-")) for x in lines]
        
        lines = [("\n  C[%d][%d][%d][%d].re = \n\t%s;\n" % (i, j, k, l, lines[l + k*NC + j*NC*NC + i*NC*NC*NS][0]),
                  "\n  C[%d][%d][%d][%d].im = \n\t%s;\n" % (i, j, k, l, lines[l + k*NC + j*NC*NC + i*NC*NC*NS][1])) 
                 for i,j,k,l in itertools.product(range(NS), range(NS), range(NC), range(NC))]
        
        buf += "".join([x[0] + x[1] for x in lines])
        buf += "\n  return;\n}\n\n"

buf += "#endif /* _QCD_PROP_CONTRACT_H */\n"
print(buf)
