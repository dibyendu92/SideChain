#-------------------------------------------------------------------------------
# . File      : ProteinResidue.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import exceptions

from ProteinContainer import ProteinContainer
from ProteinAtom import ProteinAtom
from Utilities import CalculatePositionImproper, CalculatePositionNormal


class ProteinResidue (ProteinContainer):
    def __contains__ (self, label):
        if hasattr (self, "atoms"):
            for atom in self.atoms:
                if (label == atom.label):
                    return True
        return False

    def __getitem__ (self, label):
        if hasattr (self, "atoms"):
            for atom in self.atoms:
                if (label == atom.label):
                    return atom
        raise exceptions.StandardError ("Atom %s not found in residue %s" % (label, self.Label ()))

    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0

    def _AddAtom (self, atom, build=True):
        if (build):
            if not hasattr (self, "atoms"):
                self.atoms = []
            if (atom.label in self):
                raise exceptions.StandardError ("Atom %s already exists in residue %s" % (label, self.Label ()))
            self.atoms.append (atom)

#            residue = self.parent
#            chain = residue.parent
#            protein = chain.parent
#            protein._Write ("Added atom %s to residue %s" % (atom.label, residue.Label ()))

    def _AddAtomFromIC (self, label, serial, internal, build=True):
        if (build):
            (a, b, c, distance, angle, torsion, improper) = internal
            pa = self[a].coordinates
            pb = self[b].coordinates
            pc = self[c].coordinates
            
            CalculatePosition = CalculatePositionImproper if (improper) else CalculatePositionNormal
            pd = CalculatePosition (pa, pb, pc, distance, angle, torsion)

            proteinCoordinates = self.parent.parent.coordinates
            atom = ProteinAtom (label=label, serial=serial, x=pd[0], y=pd[1], z=pd[2], parent=self, proteinCoordinates=proteinCoordinates)
            self._AddAtom (atom)


    def Label (self):
        return "%s.%s.%d" % (self.parent.label, self.label, self.serial)

    def _GetIndex (self, label):
        if (label in self):
            for (i, atom) in enumerate (self.atoms):
                if (label == atom.label):
                    return i
        return -1

    def _LabelsToIndices (self, labels):
        indices = []
        for (i, atom) in enumerate (self.atoms):
            if (atom.label in labels):
                indices.append (i)
        return indices

if (__name__ == "__main__"): pass
