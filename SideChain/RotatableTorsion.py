#-------------------------------------------------------------------------------
# . File      : ProteinModel.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import numpy, math

from ProteinContainer import ProteinContainer
from Utilities import VectorNormalize


class RotatableTorsion (ProteinContainer):
    def __len__ (self):
        return len (self.atoms)

    def Print (self):
        labels = []
        for atom in self.atoms:
            labels.append (atom.label)
        print ("Torsion: %s-... %s" % ("-".join (labels[:3]), " ".join (labels[3:])))

    def Rotate (self, degree):
        a  = self.atoms[0].coordinates
        b  = self.atoms[1].coordinates
        c  = self.atoms[2].coordinates

        ab = b - a
        j  = VectorNormalize (b - c)
        k  = VectorNormalize (numpy.cross (ab, j))
        i  = VectorNormalize (numpy.cross (j, k))

        B = numpy.empty (shape=(3, 3))  #, dtype=numpy.float
        B[0, 0] = i[0]
        B[1, 0] = i[1]
        B[2, 0] = i[2]

        B[0, 1] = j[0]
        B[1, 1] = j[1]
        B[2, 1] = j[2]

        B[0, 2] = k[0]
        B[1, 2] = k[1]
        B[2, 2] = k[2]
        Binv = numpy.linalg.inv (B)

        theta = degree * math.pi / 180.0
        cos = math.cos (theta)
        sin = math.sin (theta)

        R = numpy.identity (3)
        R[0, 0] = cos
        R[2, 0] = -sin
        R[0, 2] = sin
        R[2, 2] = cos
        C = numpy.dot (B, numpy.dot (R, Binv))

        for atom in self.atoms[3:]:
            atom.coordinates = c + numpy.dot (C, atom.coordinates - c)


if (__name__ == "__main__"): pass
