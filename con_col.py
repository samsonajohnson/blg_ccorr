#S values from http://ads.nao.ac.jp/cgi-bin/nph-iarticle_query?1993SAAOC..15....1C&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
import numpy as np

def dwarf_vi_to_ub(x):
    temp = 0.
    c = np.array([-0.3156199,\
                       0.1929464e1,\
                       -0.9275652,\
                       -0.3492653e1,\
                       0.2993744e1,\
                       0.2513098e1,\
                       -0.3020830e1,\
                       -0.9711406e-1,\
                       0.846566,\
                       -0.2171137])
               
    shfx = 1
    shfy = 1

    x += shfx
    for ind in range(len(c)):
        temp += c[ind]*x**ind
    temp += shfy
    return temp

def giant_vi_to_ub(x):
    temp = 0.
    c = np.array([-0.2158070,\
                       0.2585109e1,\
                       -0.4540624e-1,\
                       -0.4750828e1,\
                       0.3083459e1,\
                       0.5271330e1,\
                       -0.5649927e1,\
                       -0.2388506e1,\
                       0.3415714e1,\
                       0.3779219,\
                       -0.6985604])
               
    shfx = 1
    shfy = 1

    x += shfx
    for ind in range(len(c)):
        temp += c[ind]*x**ind
    temp += shfy
    return temp

def dwarf_vi_to_bv(x):
    temp = 0.
    c = np.array([-0.6865072e-1,\
                       0.8837997,\
                       -0.3889774,\
                       -0.4998126e-2,\
                       0.3867544,\
                       -0.5422331,\
                       -0.8926476e-1,\
                       0.5194797,\
                       -0.2044681,\
                       -0.1009025,\
                       0.9543256e-1,\
                       -0.2567529e-1,\
                       0.2393742e-2])
               
    shfx = 1
    shfy = 1

    x += shfx
    for ind in range(len(c)):
        temp += c[ind]*x**ind
    temp += shfy
    return temp


def giant_vi_to_bv(x):
    temp =0.
    """
    c1 = 0.122912e-1
    c2 = 0.1261086e1
    c3 = -0.4987561
    c4 = -0.1310452e1
    c5 = 0.1273935e1
    c6 = 0.2078357e1
    c7 = -0.1592421e2
    c8 = -0.1766536e1
    c9 = 0.6099412
    c10 = 0.5644602
    """
    shfx = 1
    shfy = 1

    x += shfx
    c = np.array([0.122912e-1,\
                      0.1261086e1,\
                      -0.4987561,\
                      -0.1310452e1,\
                       0.1273935e1,\
                      0.2078357e1,\
                      -0.1592421e1,\
                      -0.1766536e1,\
                       0.6099412,\
                      0.5644602])

    for ind in range(len(c)):
        temp += c[ind]*x**ind
    temp += shfy
    return temp
