import DtAlignment.DetElement as DetElement
from ROOT import TVector3

class DriftTube(DetElement):
    #TODO check if constructor should take tube center positions and rotation instead
    def __init__(self,id,vec_bottom,vec_top):
        '''Constructor. Initializes a drift tube
        
        Parameters
        ----------
        id : int
            Unique detector ID of this drift tube
        '''
        self._ID = id
        self._bottom = vec_bottom
        self._top = vec_top