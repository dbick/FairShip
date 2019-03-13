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
    x_sum = 0
    y_sum = 0
    z_sum = 0
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
    """Parse the detector id to human readable dictionary so that a specific tube can easily be
    identified and adressed out of a bigger detector arangement
    
    Parameters
    ----------
    det_id: int
        Detector ID that is to be parsed
        
    Returns
    -------
    dict
        dictionary containing the result in human readable form
    """
    result = {}
    str_id = str(det_id)
    last_four = int(str_id[-4:]) #last four digits of the detector ID
    #parse view
    view = "X"
    result['station'] = det_id / 10000000
    if result['station'] == 1:
        if str_id[1] == '1':
            view = "U"
    elif result['station'] == 2:
        if str_id[1] == '1':
            view = "V"
            
    #parse module
    module = "T" + str(result['station'])
    if result['station'] >= 3:
        if last_four <= 12:
            module += "a"
        elif last_four <= 24:
            module += "b"
        elif last_four <= 36:
            module += "c"
        else:
            module += "d"
            
    module += view
    result['module'] = module
    result['tube'] = 0
    return result

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

def distance_to_wire(tube,mom=None,pos=None):
    """Calculates the distance of closest approach for a track and a tube.
    Note: This distance is positive if a valid track was used.
    
    Parameters
    ----------
    mom: ROOT.TVector3
        Momentum (a.k.a direction) of the track
    
    pos: ROOT.TVector3
        Position on the track
        
    Returns
    -------
    float
        Closest distance between track and wire
    """
    vtop,vbot = tube.wire_end_positions()    
    normal_vector = mom.Cross(vtop-vbot)
    vec_any_two_points = vbot - pos
    distance = abs(vec_any_two_points.Dot(normal_vector)) / normal_vector.Mag() * u.mm

    return distance
    

