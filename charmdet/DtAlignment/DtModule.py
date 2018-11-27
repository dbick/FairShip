import numpy as np
import DtAlignment.DriftTube as Tube
import DtAlignment.DetElement as DetElement

class DtModule(DetElement):
    '''Drift tube module: Physical assembly of 48 drift tubes
    Defines a module of drift tubes which resembles a physical assembly of
    drift tubes. These are arranged in four layers of twelve tubes each.
    '''
    def __init__(self,list_of_tubes,x,y,z,theta,phi):
        '''Constructor
            
        Initializes a drift tube module from the list of tubes
        that are part of this module as well as the coordinates x,y and z of
        the geometric center of this module and two angles of rotation w.r.t
        a global coordinate axis system.
        
        Parameters
        ----------
        list_of_tubes : list
            list of drift tubes that are part of this module
            
        x : np.float64
            coordinate of the module center in x direction
        
        y : np.float64
            coordinate of the module center in y direction
            
        z : np.float64
            coordinate of the module center in z direction
            
        theta : np.float64
            rotation of the module w.r.t the global z axis
            
        phi : np.float64
            rotation of the module w.r.t the global y axis
            
        '''
        super().__init__(x,y,z,theta,phi)
        self._is_aligned = False
        #TODO consider deep copy
        self._list_of_tubes = list_of_tubes
       
        
    def apply_translation(self,dx,dy,dz):
        '''Apply a translation to the whole module
        
        Applies a translation to the whole module (a.k.a for all tubes in the module)
        assuming that the relative positions of the contained drift tubes remain
        unchanged.
        
        Parameters
        ----------
        dx : np.float64
            Translation in x direction
            
        dy : np.float64
            Translation in y direction
            
        dz : np.float64
            Translation in z direction
        '''
        super().apply_translation(self,dx,dy,dz)
        for tube in self._list_of_tubes:
            tube.apply_translation(dx,dy,dz)
        
    def apply_rotation(self,dTheta,dPhi):
        '''Apply rotation to the whole module
        Applies a translation to the whole module assuming that the relative
        orientations of the individual tubes in the module remain unchanged.
        
        Parameters
        ----------
        dTheta : np.float64
            Change of rotation around the z axis
        
        dPhi : np.float64
            Change of rotation around the y axis
        '''
        
        super().apply_rotation(self,dTheta,dPhi)
        for tube in self._list_of_tubes:
            tube.apply_rotation(dTheta,dPhi)
        
    def get_tubes(self):
        '''Get a list of tubes in this module
        Returns a list of drift tubes that are part of this module
        
        Returns
        -------
        list
            List of tubes contained in this module
        '''
        return self._list_of_tubes