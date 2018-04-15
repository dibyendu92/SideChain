#-------------------------------------------------------------------------------
# . File      : ProteinContainer.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------

class ProteinContainer (object):
    def __init__ (self, **kwargs):
        for (key, value) in kwargs.items ():
            setattr (self, key, value)

if (__name__ == "__main__"): pass
