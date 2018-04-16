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
        return self._coordinates[self.serial-1, :]

    @coordinates.setter
    def coordinates (self, array):
        for i in range (3):
            self._coordinates[self.serial-1, i] = array[i]

    @property
    def x (self):
        return self._coordinates[self.serial-1, 0]

    @x.setter
    def x (self, component):
        self._coordinates[self.serial-1, 0] = component

    @property
    def y (self):
        return self._coordinates[self.serial-1, 1]

    @y.setter
    def y (self, component):
        self._coordinates[self.serial-1, 1] = component

    @property
    def z (self):
        return self._coordinates[self.serial-1, 2]

    @z.setter
    def z (self, component):
        self._coordinates[self.serial-1, 2] = component

    def __init__ (self, label, serial, x, y, z, parent, proteinCoordinates):
        self.label = label
        self.serial = serial
        self.parent = parent
        self._coordinates = proteinCoordinates

        if (self._coordinates != None):
            component = (x, y, z)
            for i in range (3):
                self._coordinates[self.serial-1, i] = component[i]

if (__name__ == "__main__"): pass
