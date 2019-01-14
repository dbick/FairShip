from ROOT import TVector3, TRotation
import shipunit as u

def calculate_center_from_lot(list_of_tubes):
    return -1
    
def calculate_rotation_from_lot(list_of_tubes):
    return -1 
    
def parse_det_id(det_id):
    return -1

def calculate_center(vec1, vec2):
    '''Calculate the position in the middle between two positions pointed to by vectors.
    
    Parameters
    ----------
    vec1: TVector3
        Vector to the first position
        
    vec2: TVector3
        Vector to the second position
        
    Returns
    -------
    TVector3
        Position of the point in between the two positions passed as arguments
    '''
    vec_between = vec2 - vec1
    half_length = vec_between.Mag() / 2
    vec_between.SetMag(half_length)
    
    return TVector3(vec1 + vec_between)
    
def z_rotation_to_euler_angles(rad_z):
    '''Convert a rotation around the z-axis to Euler angles following the X-convention.
    
    Parameters
    ----------
    rad_z: float
        Radians of rotation around z
        
    Returns
    -------
    float,float,float
        Phi, Theta, Psi following the Euler-X-convention
    '''
    rot = TRotation()
    rot.RotateZ(rad_z)    
    return rot.GetXPhi(), rot.GetXTheta(), rot.GetXPsi()
