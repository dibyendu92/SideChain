#-------------------------------------------------------------------------------
# . File      : RotatableBonds.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, exceptions


class RotatableBonds (object):
    """A class collecting rotatable bonds for every residue."""

    def __init__ (self, filename=os.path.join ("toppar", "rotatable.dat"), verbose=True):
        """Constructor."""
        self.prompt = "# . RotatableBonds> "
        self.verbose = verbose
        self.pairs = {}

        lines = open (filename).readlines ()
        for line in lines:
            tokens = line.split ()
            component = tokens[0]
            labels = tokens[1:]
            nlabels = len (labels)
            if ((nlabels % 2) != 0):
                raise exceptions.StandardError ("Incorrect rotatable bond format for component %s." % component)
            pairs = []
            for i in range (0, nlabels, 2):
                pair = (labels[i], labels[i + 1])
                pairs.append (pair)
            self.pairs[component] = pairs
        self._Write ("Read %s components." % len (self))

    def _Write (self, message):
        if (self.verbose):
            if hasattr (self, "prompt"):
                print ("%s%s" % (self.prompt, message))

    def __contains__ (self, component):
        return self.pairs.has_key (component)

    def __getitem__ (self, component):
        if (component in self):
            return self.pairs[component]
        raise exceptions.StandardError ("Rotatable bonds for comonent %s not found." % component)

    def __len__ (self):
        return len (self.pairs)

#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
