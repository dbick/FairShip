import numpy as np

#TODO: Two rotational angles only applicable to vectors - consider using Euler angles for Module

class Abstract_DetElement:
    '''Detector element abstract class. Do not instantiate.
    Contains basic members and methods needed for every detector element
    '''
    
    def __init__(self,x,y,z,theta,phi,rho):
        '''Constructor
            
        Initializes a detector element from the coordinates x,y and z of
        the geometric center of this element and two angles of rotation w.r.t
        a global coordinate axis system.
        
        Parameters
        ----------
        x : np.float64
            coordinate of the element's center in x direction
        
        y : np.float64
            coordinate of the element's center in y direction
            
        z : np.float64
            coordinate of the element's center in z direction
            
        theta : np.float64
            rotation of the element w.r.t the global z axis
            
        phi : np.float64
            rotation of the element w.r.t the global y axis
        
        rho : np.float64
            rotation of the element w.r.t the global x axis
        '''
        self._position = np.array([x,y,z], dtype=np.float64)
        self._rotation = np.array([theta,phi,rho], dtype = np.float64)
                
    def apply_translation(self,dx,dy,dz):
        '''Apply a translation to the whole module
        
        Applies a translation to the element
                
        Parameters
        ----------
        dx : np.float64
            Translation in x direction
            
        dy : np.float64
            Translation in y direction
            
        dz : np.float64
            Translation in z direction
        '''
        self._position[0] += dx
        self._position[1] += dy
        self._position[2] += dz
        
    def apply_rotation(self,dTheta,dPhi,dRho):
        '''Apply rotation to the detector element
        Applies a translation to the element.
        
        Parameters
        ----------
        dTheta : np.float64
            Change of rotation around the z axis
        
        dPhi : np.float64
            Change of rotation around the y axis
        '''
        
        self._rotation[0] += dTheta
        self._rotation[1] += dPhi
        self._rotation[2] += dRho