#!/usr/bin/python
import os, sys, getopt, re, struct
import itertools, time
import scipy, scipy.fftpack

default_mom_list = "mom-list"
default_suppress_idx = False

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "f:nh", ["mom-file=", "no-idx", "help"])
        except getopt.GetoptError as msg:
            raise Usage(msg)

        for o, a in opts:
            if o in ["-h", "--help"]:
                print(" Usage: %s [OPTIONS] FILE" % argv[0], file=sys.stderr)
                print(" Options:", file=sys.stderr)
                print("  -f, --mom-file\t\t: Momentum list file", file=sys.stderr)
                print("  -n, --no-idx\t\t: Suppress output of last index (e.g. for scalar, pseudo-scalar where there is only a single entry", file=sys.stderr)
                print("  -h, --help\t\t: This help message", file=sys.stderr)
                print 
                return 2

            elif o in ["-f", "--mom-file"]:
                user_mom_list = a

            elif o in ["-n", "--no-idx"]:
                user_suppress_idx = True

            else:
                print(" %s: ignoring unhandled option" % o, file=sys.stderr)

        if len(args) != 1:
            raise Usage(" Usage: %s [OPTIONS] FILE" % argv[0])

    except Usage as err:
        print(err.msg, file=sys.stderr)
        print(" for help use --help", file=sys.stderr)
        return 2

    try:
        mom_list = user_mom_list
    except NameError:
        mom_list = default_mom_list

    try:
        suppress_idx = user_suppress_idx
    except NameError:
        suppress_idx = default_suppress_idx

    ar = scipy.array
    fname = args[0]

    fp = open(fname, "rb")
    line = fp.readline().decode()
    ### first line should contain "begin-header"
    if line != "begin-header;\n":
        print("Malformed header")
        raise IOError

    ### Now get the rest of the header
    header = [line.split(";")[0],]
    while True:
        line = fp.readline()
        if len(line) == 0:
            print("File ended unexpectedly")
            raise IOError

        line = line.decode().split(";")[0]
        header.append(line)
        if line == "end-header":
            break

    ### The data
    buf = fp.read()
    fp.close()

    ### Header parser
    def parse(tag, header):
        line = [x for x in header if re.search(tag, x)]
        if len(line) != 1:
            print("%s not found in header" % tag)
            raise ValueError
    
        return line[0].split("=")[1].strip()


    dims = [int(x) for x in parse("lattice dimensions", header).split(",")]
    src_pos = [int(x) for x in parse("source position", header).split(",")]
    site_tags = [x.strip() for x in parse("lattice-site tags", header).split(",")]
    site_size = len(site_tags)
    t_ins = [int(x) for x in parse("t-insertions", header).split(",")]
    nt_ins = len(t_ins)

    nelems = site_size*nt_ins*dims[0]*dims[1]*dims[2]*2
    if len(buf) != nelems*8:
        print("binary data size miss-match")
        raise ValueError

    fmt = ">" + "".join(nelems*["d"])
    buf = struct.unpack(fmt, buf)
    buf = ar(buf).reshape([-1, 2])
    thrp = buf[:, 0] + complex(0, 1)*buf[:, 1]
    thrp = thrp.reshape([nt_ins, dims[2], dims[1], dims[0], site_size])

    for ax,shift in zip([0, 1, 2], [src_pos[2],src_pos[1],src_pos[0]]):
        thrp = scipy.roll(thrp, shift=-shift, axis=ax+1)

    thrp_p = []
    for t in range(nt_ins):
        p = []
        for s in range(site_size):
            x = thrp[t, :, :, :, s]
            p.append(scipy.fftpack.fftn(x))
        thrp_p.append(p)

    moms = open(mom_list, "r").readlines()
    moms = [x.split() for x in moms]
    moms = [[int(x) for x in y] for y in moms]
    for it,t in enumerate(t_ins):
        for k,s in enumerate(site_tags):
            name = s.replace("(", "").replace(")", "").split(":")[0]
            index = s.replace("(", "").replace(")", "").split(":")[1].strip()
            for m in moms:
                x = (dims[0] - m[0]) % dims[0]
                y = (dims[1] - m[1]) % dims[1]
                z = (dims[2] - m[2]) % dims[2]
                if suppress_idx:
                    print("%s %+d %+d %+d %+13.6e %+13.6e" % (t, m[0], m[1], m[2], 
                                                              thrp_p[it][k][z, y, x].real, 
                                                              thrp_p[it][k][z, y, x].imag))
                else:
                    print("%s %+d %+d %+d %+13.6e %+13.6e %s" % (t, m[0], m[1], m[2], 
                                                                 thrp_p[it][k][z, y, x].real, 
                                                                 thrp_p[it][k][z, y, x].imag, index))
                    
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
