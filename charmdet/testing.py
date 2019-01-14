#Tested in Python 2.7.6
import numpy as np
import DtAlignment.DriftTube as Tube
import DtAlignment.DtModule as DtModule
import DtAlignment.utils
from ROOT import TH3D
from ROOT import TVector3
import shipunit as u

#create 1m tubes in 4 layers of twelve:
#assume all tube centers are at y = 0
test_id = 10000000
list_of_tubes = []

nBins = 1000000
bin_range = [-1*u.m, 1*u.m]
center_positions = TH3D("Positions","Tube center positions",nBins,bin_range[0],bin_range[1],nBins,bin_range[0],bin_range[1],nBins,bin_range[0],bin_range[1])
x = 0.0
y = 0
z_layer = 0
for layer in range(4):
    z_layer = layer * 4.2
    for i in range(12):
        x += 4.2
        tube_id = (test_id + i) + (10 * layer)
        list_of_tubes.append(Tube(tube_id,x,y,z_layer,0,0,0))
    x = 0
    y = 0

module = DtModule(list_of_tubes,12 * 4.2 / 2,0,0,0,0,0)
#remember data at 0 rotatation
tube_xyz_rot0 = []
for tube in module.get_tubes():
    xyz = [np.float64(tube._position[0]),np.float64(tube._position[1]),np.float64(tube._position[2])]
    #TODO append a COPY of xyz to tube_xyz_rot0
    tube_xyz_rot0.append(xyz)
    
    
print("Table: x,y,z for each tube")
for tube in module.get_tubes():
    x = tube._position[0]
    y = tube._position[1]
    z = tube._position[2]
    print("{}\t{}\t{}".format(x,y,z))
    #center_positions.Fill(np.float64(tube._position[0]),np.float64(tube._position[1]),np.float64(tube._position[2]))
    
print("Rotated")
module.rotate_z(60 * 2 * np.pi / 360)
for tube in module.get_tubes():
    x = tube._position[0]
    y = tube._position[1]
    z = tube._position[2]
    print("{}\t{}\t{}".format(x,y,z))
    
print("Translated")
module.apply_translation(-5*u.cm,2 * u.cm, 3.7*u.mm)
for tube in module.get_tubes():
    x = tube._position[0]
    y = tube._position[1]
    z = tube._position[2]
    print("{}\t{}\t{}".format(x,y,z))

print(DtAlignment.utils.calculate_center_from_lot([]))
    