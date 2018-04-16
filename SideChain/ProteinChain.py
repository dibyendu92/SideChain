#-------------------------------------------------------------------------------
# . File      : ProteinChain.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
from ProteinContainer import ProteinContainer


class ProteinChain (ProteinContainer):
    def _AddResidue (self, residue, build=True):
        if (build):
            if not hasattr (self, "residues"):
                self.residues = []
            self.residues.append (residue)

if (__name__ == "__main__"): pass
