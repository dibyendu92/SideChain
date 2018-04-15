#-------------------------------------------------------------------------------
# . File      : InternalLibrary.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import os, glob, exceptions

from MolarisTools.Library import InternalCoordinate
from MolarisTools.Utilities import TokenizeLine


class InternalLibrary (object):
    """A class collecting tables of internal coordinates."""

    @staticmethod
    def ReadICFile (filename):
        lines = open (filename, "r").readlines ()
        ictable = []
        (label, natoms) = TokenizeLine (lines[0], converters=[None, int])
        for line in lines[1:]:
            tokens = TokenizeLine (line, converters=[None, ] * 4 + [float, ] * 5 + [int, ])
            ic = InternalCoordinate (i = tokens[0], 
                                    j = tokens[1], 
                                    k = tokens[2], 
                                    l = tokens[3], 
                                    distanceLeft = tokens[4], 
                                    angleLeft = tokens[5], 
                                    torsion = tokens[6], 
                                    angleRight = tokens[7], 
                                    distanceRight = tokens[8], 
                                    improper = (True if (tokens[9] == 1) else False) )
            ictable.append (ic)    
        return (label, ictable)

    def __init__ (self, path="toppar", verbose=True):
        """Constructor."""
        self.path = path
        self.verbose = verbose
        self.prompt = "# . InternalLibrary> "
        self._ReadFiles ()

    def _ReadFiles (self):
        self.components = {}
        files = glob.glob (os.path.join (self.path, "*.ic"))
        for filename in files:
            (label, ictable) = self.ReadICFile (filename)
            self.components[label] = ictable
            self._Write ("Read file \"%s\"" % filename)

    def _Write (self, message):
        if (self.verbose):
            if hasattr (self, "prompt"):
                print ("%s%s" % (self.prompt, message))

    def __len__ (self):
        if hasattr (self, "components"):
            return len (self.components)
        return 0

    def __contains__ (self, key):
        if hasattr (self, "components"):
            return self.components.has_key (key)
        return False

    def __getitem__ (self, key):
        if hasattr (self, "components"):
            if (key in self):
                return self.components[key]
        raise exceptions.StandardError ("Internal coordinates for component %s not found." % key)


if (__name__ == "__main__"): pass
