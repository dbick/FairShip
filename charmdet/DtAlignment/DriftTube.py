import DtAlignment.DetElement as DetElement
from ROOT import TVector3
import shipunit as u

class DriftTube(DetElement):
    #TODO check if constructor should take tube center positions and rotation instead
    def __init__(self,detId,vec_bottom,vec_top):
        '''Constructor. Initializes a drift tube.
        
        Parameters
        ----------
        id : int
            Unique detector ID of this drift tube
        '''
        #TODO: Read length from TGeo
        if(int(detId / 1000000) < 3):
            self.length = 100 * u.cm
        else:
            self._length = 160 * u.cm
            
        self._ID = detId
        self._bottom = vec_bottom
        self._top = vec_top
        
    def wire_end_positions(self):
        '''Returns the end positions of the wire in a global reference system.
        
        Returns
        -------
        TVector3
            Top position of the wire
            
        TVector3
            Bottom position of the wire
        '''
        l = self._length / 2
        x_top = self._position[0] #TODO consider rotation
        y_top = self._position[1] + l #TODO consider rotation
        z_top = self._position[2] #TODO consider rotation
        
        x_bot = self._position[0] #TODO consider rotation
        y_bot = self._position[1] - l #TODO consider rotation
        z_bot = self._position[2] #TODO consider rotation
        
        return TVector3(x_top,y_top,z_top), TVector3(x_bot,y_bot,z_bot)