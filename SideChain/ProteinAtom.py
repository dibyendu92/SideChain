#-------------------------------------------------------------------------------
# . File      : ProteinAtom.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import numpy

from ProteinContainer import ProteinContainer


class ProteinAtom (ProteinContainer):
    @property
    def coordinates (self):
        return numpy.array ((self.x, self.y, self.z))

    @coordinates.setter
    def coordinates (self, array):
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]

if (__name__ == "__main__"): pass
