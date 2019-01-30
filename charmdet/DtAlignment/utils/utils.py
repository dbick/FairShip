from ROOT import TVector3, TRotation
import shipunit as u

def calculate_center_from_lot(list_of_tubes):
    '''Calculates the geometric center point of a DtModule from a given list of DriftTube objects
    that are passed to this function as parameter list_of_tubes.
    
    Parameters
    ----------
    list_of_tubes: list
        List of drift tube objects that are part of the module for which the geometric center is to be computed
        
    Returns
    -------
    TVector3
        Vector to the center point of the module
    '''
    x_sum, y_sum, z_sum = 0,0,0
    for tube in list_of_tubes:
        tube_center = tube.get_center_position()
        x_sum += tube_center[0]
        y_sum += tube_center[1]
        z_sum += tube_center[2]
    
    center_x = x_sum / len(list_of_tubes)
    center_y = y_sum / len(list_of_tubes)
    center_z = z_sum / len(list_of_tubes)
    
    return TVector3(center_x,center_y,center_z)
    
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
