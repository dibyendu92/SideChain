#-------------------------------------------------------------------------------
# . File      : Utilities.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import numpy, math

def VectorNormalize (vector):
    return (1.0 / numpy.linalg.norm (vector)) * vector

def CalculatePositionImproper (a, b, c, Rcd, Tbcd, Pabcd):
#                I (a)    L (d)
#                 \     /
#                  \   /
#                   *K (c)
#                   |
#                   |
#                   J (b)
#        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)
    ca = a - c
    j  = VectorNormalize (b - c)
    k  = VectorNormalize (numpy.cross (ca, j))
    i  = VectorNormalize (numpy.cross (j, k))

    psi = Pabcd * math.pi / 180.0
    t   = math.cos (psi) * i + math.sin (psi) * k

    chi = Tbcd * math.pi / 180.0
    q   = math.cos (chi) * j + math.sin (chi) * t

    return (c + Rcd * q)

def CalculatePositionNormal (a, b, c, Rcd, Tbcd, Pabcd):
#            (a) I
#                 \
#                  \
#               (b) J----K (c)
#                         \
#                          \
#                           L (d)
#        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)
    ba  = a - b
    j   = VectorNormalize (c - b)
    k   = VectorNormalize (numpy.cross (ba, j))
    i   = VectorNormalize (numpy.cross (j, k))

    psi = Pabcd * math.pi / 180.0
    t   = math.cos (psi) * i + math.sin (psi) * k

    chi = Tbcd * math.pi / 180.0
    q   = -math.cos (chi) * j + math.sin (chi) * t

    return (c + Rcd * q)


if (__name__ == "__main__"): pass
