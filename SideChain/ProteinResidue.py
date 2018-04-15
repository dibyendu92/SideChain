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
        raise exceptions.StandardError ("Atom %s not found in residue %s.%d.%s" % (label, self.parent.label, self.serial, self.label))

    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0

    def AddAtom (self, atom):
        if not hasattr (self, "atoms"):
            self.atoms = []
        if (atom.label not in self):
            self.atoms.append (atom)

    def AddAtomFromIC (self, label, serial, internal):
        (a, b, c, distance, angle, torsion, improper) = internal

        pa = self[a].coordinates
        pb = self[b].coordinates
        pc = self[c].coordinates
        
        CalculatePosition = CalculatePositionImproper if (improper) else CalculatePositionNormal
        pd = CalculatePosition (pa, pb, pc, distance, angle, torsion)
        
        atom = ProteinAtom (label=label, serial=serial, x=pd[0], y=pd[1], z=pd[2], parent=self)
        self.AddAtom (atom)


    def Label (self):
        return "%s.%s.%d" % (self.parent.label, self.label, self.serial)

    def GetIndex (self, label):
        if (label in self):
            for (i, atom) in enumerate (self.atoms):
                if (label == atom.label):
                    return i
        return -1

    def LabelsToIndices (self, labels):
        indices = []
        for (i, atom) in enumerate (self.atoms):
            if (atom.label in labels):
                indices.append (i)
        return indices

if (__name__ == "__main__"): pass
