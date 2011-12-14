#!/usr/bin/python
import os, sys, struct
import scipy, pylab

### some local defs
dot = scipy.dot
sum = scipy.sum
sqrt = scipy.sqrt

### Lattice length
Ls = 48
fname = ["block_a5p0_N220.dat", "block_a0p1_N420.dat"]

dat = []
for f in fname:
### open and read binary file, store as array
    buf = open(f, "br").read()
    Vol = Ls**3
    fmt = ">" + "".join(["d"]*Vol)
    s = struct.Struct(fmt)
    buf = struct.unpack(fmt, buf)
    buf = scipy.array(buf).reshape([Ls, Ls, Ls])
    zyx = scipy.mgrid[:Ls,:Ls,:Ls].reshape(3, -1)
    dat.append(buf)
### loop over sites to compute r.m.s radius
for d in dat:
    rms = 0
    norm = 0
    for i in zip(*zyx):
        iz, iy, ix = i
        x = ix if ix < Ls/2 else Ls - ix
        y = iy if iy < Ls/2 else Ls - iy
        z = iz if iz < Ls/2 else Ls - iz
        rms += (x*x + y*y + z*z) * d[iz, iy, ix]
        norm += d[iz, iy, ix]
    print(" r^2 = %e" % (rms/norm))

### loop over sites to compute profile 
func_r = []
count_r = []
for d in dat:
    f_r = {}
    c_r = {}
    for i in zip(*zyx):
        iz, iy, ix = i
        x = ix if ix < Ls/2 else Ls - ix
        y = iy if iy < Ls/2 else Ls - iy
        z = iz if iz < Ls/2 else Ls - iz
        r2 = x*x + y*y + z*z
        psi = d[iz, iy, ix] 
        try:
            f_r[r2] += psi
            c_r[r2] += 1
        except KeyError:
            f_r[r2] = psi
            c_r[r2] = 1
        norm += psi

### normalize profile
    for r2 in f_r:
        f_r[r2] /= (norm*c_r[r2])

    func_r.append(f_r)
    count_r.append(c_r)

fig = pylab.figure(1)
fig.clf()
ax = fig.gca()
ax.set_yscale("log")
symb = ["ro", "bs"]
for s, name, f in zip(symb, fname, func_r):
### plot profile
    xy = [(r2, f[r2]) for r2 in f]
    xy = [x for x in zip(*sorted(xy))]
    ax.plot(sqrt(xy[0]), xy[1], s, label=name)
ax.legend()
fig.canvas.draw()
fig.show()
