import DtAlignment.DriftTube as DriftTube
import DtAlignment.DtModule as DtModule
import DtAlignment.utils
import ROOT
import shipunit as u

xpos = {}
xposb = {}
ypos = {}
yposb = {}
zpos = {}
residuals = [0.]*24
# positions are relative to the top / bottom end plates of a station, corrected from survey positions with known offset in y, 
# 2cm + length of the bolt 150mm on the top and 50mm on the bottom
Nchannels = {1:12,2:12,3:48,4:48}

#survey
survey = {} # X survey[xxx][1] Y survey[xxx][2] Z survey[xxx][0]
daniel = {}
# For T1 and 2, the survey target is placed above/below the endplates, offset in y. 
# For T3 it is placed 7cm in front and for T4 7cm behind the endplate, so there is an offset in z
survey['T1_MA_01']=[ 9.0527, 0.2443, 0.7102 ]
daniel['T1_MA_01']=[244.30,531.85,9052.70]
survey['T1_MA_02']=[ 9.0502, -0.2078, 0.7092  ]
daniel['T1_MA_02']=[-207.80,530.85,9050.20]
survey['T1_MA_03']=[ 9.0564, -0.2075, -0.6495  ]
daniel['T1_MA_03']=[-207.50,-570.85,9056.40]
survey['T1_MA_04']=[ 9.0578, 0.2436, -0.6503  ]
daniel['T1_MA_04']=[ 243.60,-571.65,9057.80]
survey['T1_MB_01']=[ 9.5376, -0.5044, 0.5647  ]
daniel['T1_MB_01']=[ -504.54,564.60,9537.42  ]
survey['T1_MB_02']=[ 9.5388, -0.7293, 0.1728  ]
daniel['T1_MB_02']=[ -729.38,172.83,9539.00  ]
survey['T1_MB_03']=[ 9.5420, 0.4446, -0.4995  ]
daniel['T1_MB_03']=[ 444.63,-499.49,9541.94  ]
survey['T1_MB_04']=[ 9.5435, 0.6693, -0.1073  ]
daniel['T1_MB_04']=[ 669.28,-107.28,9543.47  ]
survey['T2_MC_01']=[ 9.7424, 0.7588, 0.1730  ]
daniel['T2_MC_01']=[ 758.84, 173.03, 9742.35  ]
survey['T2_MC_02']=[ 9.7433, 0.5327, 0.5657  ]
daniel['T2_MC_02']=[532.67, 565.72, 9743.27 ]
survey['T2_MC_03']=[ 9.7381, -0.6374, -0.1104  ]
daniel['T2_MC_03']=[-637.45, -110.36, 9738.21 ]
survey['T2_MC_04']=[ 9.7409, -0.4121, -0.5021  ]
daniel['T2_MC_04']=[-412.30,-501.89, 9740.98  ]
survey['T2_MD_01']=[ 10.2314, 0.2385, 0.7078  ]
daniel['T2_MD_01']=[ 238.50, 530.50, 10231.40  ]
survey['T2_MD_02']=[ 10.2276, -0.2121, 0.7087  ]
daniel['T2_MD_02']=[ -212.10,531.40,10227  ]
survey['T2_MD_03']=[ 10.2285, -0.2157, -0.6488  ]
daniel['T2_MD_03']=[ -215.70, -570.35, 10228.50 ]
survey['T2_MD_04']=[ 10.2302, 0.2361, -0.6495  ]
daniel['T2_MD_04']=[ 236.10, -575.05, 10230.20  ]
survey['T3_B01']=[  14.5712,  0.9285, -0.6818 ]
daniel['T3_B01']=[  928.59, -681.74, 14641.10 ]
survey['T3_B02']=[  14.5704,  0.5926, -0.6823 ]
daniel['T3_B02']=[  592.62, -682.23, 14640.30 ]
survey['T3_B03']=[  14.5699,  0.4245, -0.6844 ]
daniel['T3_B03']=[  424.52, -684.35, 14639.90 ]
survey['T3_B04']=[  14.5686,  0.0884, -0.6854 ]
daniel['T3_B04']=[  88.46, -685.39, 14638.60 ]
survey['T3_B05']=[  14.5685, -0.0813, -0.6836 ]
daniel['T3_B05']=[  -81.23, -683.57, 14638.50 ]
survey['T3_B06']=[  14.5694, -0.4172, -0.6840 ]
daniel['T3_B06']=[  -417.09, -683.99, 14639.30 ]
survey['T3_B07']=[  14.5696, -0.5859, -0.6864 ]
daniel['T3_B07']=[  -585.80, -686.33, 14639.60 ]
survey['T3_B08']=[  14.5693, -0.9216, -0.6845 ]
daniel['T3_B08']=[  -921.58, -684.47, 14639.20 ]
survey['T3_T01']=[  14.5733,  0.9253, 0.8931 ]
daniel['T3_T01']=[  925.40, 893.09,14643.20 ]
survey['T3_T02']=[  14.5741,  0.5893, 0.8914 ]
daniel['T3_T02']=[  589.42, 891.36,14644.10 ]
survey['T3_T03']=[  14.5746,  0.4212, 0.8907 ]
daniel['T3_T03']=[  421.23, 890.75, 14644.60 ]
survey['T3_T04']=[  14.5750,  0.0852, 0.8905 ]
daniel['T3_T04']=[  85.28, 890.55, 14645.00 ]
survey['T3_T05']=[  14.5756, -0.0839, 0.8899 ]
daniel['T3_T05']=[  -83.89, 889.94, 14645.50 ]
survey['T3_T06']=[  14.5769, -0.4198, 0.8888 ]
daniel['T3_T06']=[  -419.69, 888.85, 14646.90 ]
survey['T3_T07']=[  14.5781, -0.5896, 0.8908 ]
daniel['T3_T07']=[  -589.56, 890.77, 14648.10 ]
survey['T3_T08']=[  14.5812, -0.9256, 0.8896 ]
daniel['T3_T08']=[  -925.57, 889.62, 14651.20 ]
survey['T4_B01']=[  16.5436,  0.9184, -0.6848 ]
daniel['T4_B01']=[  918.35,-684.86,16473.60 ]
survey['T4_B02']=[  16.5418,  0.5824, -0.6867 ]
daniel['T4_B02']=[  582.36, -686.73, 16471.90 ]
survey['T4_B03']=[  16.5408,  0.4144, -0.6875 ]
daniel['T4_B03']=[  414.37, -687.52, 16470.90 ]
survey['T4_B04']=[  16.5389,  0.0785, -0.6883 ]
daniel['T4_B04']=[  78.51, -688.32, 16468.90 ]
survey['T4_B05']=[  16.5389, -0.0888, -0.6890 ]
daniel['T4_B05']=[  -88.80, -689.05, 16468.90 ]
survey['T4_B06']=[  16.5396, -0.4247, -0.6884 ]
daniel['T4_B06']=[  -424.75, -688.49, 16469.60 ]
survey['T4_B07']=[  16.5400, -0.5936, -0.6888 ]
daniel['T4_B07']=[  -593.65, -688.85, 16470.00 ]
survey['T4_B08']=[  16.5402, -0.9295, -0.6877 ]
daniel['T4_B08']=[  -929.54, -687.77, 16470.20 ]
survey['T4_T01']=[  16.5449,  0.9207, 0.8899 ]
daniel['T4_T01']=[  920.70,889.81,16475.00  ]
survey['T4_T02']=[  16.5456,  0.5845, 0.8884 ]
daniel['T4_T02']=[  584.44, 888.30, 16475.70  ]
survey['T4_T03']=[  16.5460,  0.4168, 0.8873 ]
daniel['T4_T03']=[  416.74, 887.26, 16476.10  ]
survey['T4_T04']=[  16.5474,  0.0804, 0.8862 ]
daniel['T4_T04']=[  80.34, 886.14, 16477.50 ]
survey['T4_T05']=[  16.5480, -0.0876, 0.8862 ]
daniel['T4_T05']=[  -87.55, 886.10, 16478.00 ]
survey['T4_T06']=[  16.5497, -0.4234, 0.8862 ]
daniel['T4_T06']=[  -423.40, 886.18, 16479.80 ]
survey['T4_T07']=[  16.5507, -0.5919, 0.8866 ]
daniel['T4_T07']=[  -591.89, 886.58, 16480.80 ]
survey['T4_T08']=[  16.5518, -0.9280, 0.8869 ]
daniel['T4_T08']=[  -928.06, 886.85, 16481.80 ]

survey['RPC1_L']= [ 17.6823, 1.1611, 1.1909]
survey['RPC1_R']= [ 17.6864, -1.2679, 1.2145 ]
survey['RPC2_L']= [ 18.6319, 1.1640, 1.1926 ]
survey['RPC2_R']= [ 18.6360, -1.2650, 1.2065 ]
survey['RPC3_L']= [ 19.1856, 1.1644, 1.1933 ]
survey['RPC3_R']= [ 19.1902, -1.2646, 1.2021 ]
survey['RPC4_L']= [ 19.7371, 1.1610, 1.1938 ]
survey['RPC4_R']= [ 19.7410, -1.2670, 1.1979 ]
survey['RPC5_L']= [ 20.2852, 1.1677, 1.1945 ]
survey['RPC5_R']= [ 20.2891, -1.2614, 1.1943 ]

Lcorrection={}
# length of the bolt 150mm on the top and 50mm on the bottom
Lcorrection['T1_MA_01'] = -(20.+150.+8.65)/10.
Lcorrection['T1_MA_02'] = -(20.+150.+8.65)/10.
Lcorrection['T1_MA_03'] = (20.+50.+8.35)/10.
Lcorrection['T1_MA_04'] = (20.+50.+8.35)/10.
Lcorrection['T1_MB_01'] = -(20.+150.)/10.
Lcorrection['T1_MB_02'] = -(20.+150.)/10.
Lcorrection['T1_MB_03'] = (20.+50.)/10.
Lcorrection['T1_MB_04'] = (20.+50.)/10.
Lcorrection['T2_MC_01'] = -(20.+150.)/10.
Lcorrection['T2_MC_02'] = -(20.+150.)/10.
Lcorrection['T2_MC_03'] = (20.+50.)/10.
Lcorrection['T2_MC_04'] = (20.+50.)/10.
Lcorrection['T2_MD_01'] = -(20.+150.+7.3)/10.
Lcorrection['T2_MD_02'] = -(20.+150.+7.3)/10.
Lcorrection['T2_MD_03'] = (20.+50.+8.45)/10.   # 8.45 or 4.45
Lcorrection['T2_MD_04'] = (20.+50.+8.45)/10.   # 8.45 or 4.45

surveyXYZ = {}
for s in survey: 
  surveyXYZ[s]=[survey[s][1]*100.,survey[s][2]*100.,survey[s][0]*100.]
for s in daniel: daniel[s]=[daniel[s][0]/10.,daniel[s][1]/10.,daniel[s][2]/10.]

Langle={}
for s in ['T1_MA_','T1_MB_','T2_MC_','T2_MD_']:
   delx = surveyXYZ[s+'01'][0]-surveyXYZ[s+'04'][0]
   dely = surveyXYZ[s+'01'][1]-surveyXYZ[s+'04'][1]
   Langle[s+'01-04'] = ROOT.TMath.ATan2(dely,delx)
   delx = surveyXYZ[s+'02'][0]-surveyXYZ[s+'03'][0]
   dely = surveyXYZ[s+'02'][1]-surveyXYZ[s+'03'][1]
   Langle[s+'02-03'] = ROOT.TMath.ATan2(dely,delx)
for s in survey: 
  if s.find('T1')<0 and s.find('T2')<0: continue
  if s.find('01')>0:   a = s.replace('01','01-04')
  elif s.find('04')>0: a = s.replace('04','01-04')
  elif s.find('02')>0: a = s.replace('02','02-03')
  elif s.find('03')>0: a = s.replace('03','02-03')
  #surveyXYZ[s][0] = surveyXYZ[s][0] + Lcorrection[s]*ROOT.TMath.Cos(Langle[a])
  #surveyXYZ[s][1] = surveyXYZ[s][1] + Lcorrection[s]*ROOT.TMath.Sin(Langle[a])


rn ={}
rn['T1_MB_01'] = [21.55,72.60]    # top left = 1
rn['T1_MB_02'] = [-23.65,72.60]   # top right = 2
rn['T1_MB_04'] = [21.55,-62.60]   # bottom left = 4
rn['T1_MB_03'] = [-23.65,-62.60]  # bottom right = 3
rn['T2_MC_01'] = rn['T1_MB_01'] 
rn['T2_MC_02'] = rn['T1_MB_02']
rn['T2_MC_04'] = rn['T1_MB_04']
rn['T2_MC_03'] = rn['T1_MB_03']

#Stefan: Build set of DT modules
dt_modules = {}
tubes = {}

tubes['T1X'] = []
#overall z positioning
#T1X:
zpos['T1X'] = (daniel['T1_MA_01'][2]+daniel['T1_MA_02'][2]+daniel['T1_MA_03'][2]+daniel['T1_MA_04'][2])/4. + 3.03
ypos['T1X'] = [(daniel['T1_MA_01'][1]+daniel['T1_MA_02'][1])/2.,(daniel['T1_MA_04'][1]+daniel['T1_MA_03'][1])/2.]

deltaZ = zpos['T1X'] 
n = 10002012
start = daniel['T1_MA_01'][0] # (daniel['T1_MA_01'][0]+daniel['T1_MA_04'][0])/2. bottom survey measurements do not match nominal distance
delta = (start - (daniel['T1_MA_02'][0]+daniel['T1_MA_03'][0])/2.-1.1-2.1) / 10.
delta = 4.2
for i in range(12): 
 xpos[n-i] = start - delta * i
 ypos[n-i] = ypos['T1X']
 zpos[n-i] = zpos['T1X']-deltaZ
 y_center = ypos[n-i][1] + ((ypos[n-i][0]-ypos[n-i][1]) / 2)
 tubes['T1X'].append(DriftTube(n-i,xpos[n-i],y_center,zpos[n-i]))
n = 10012001
start =  daniel['T1_MA_02'][0] +1.1 #   (daniel['T1_MA_02'][0]+daniel['T1_MA_03'][0])/2. +1.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T1X']
 zpos[n+i] = zpos['T1X']+3.64-deltaZ
 print(n+i,ypos[n+i])
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T1X'].append(DriftTube(n-i,xpos[n+i],y_center,zpos[n+i]))
n = 10102001
start = start -1.1 #  (daniel['T1_MA_02'][0]+daniel['T1_MA_03'][0])/2.
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T1X']
 zpos[n+i] = zpos['T1X']+3.64+4.06-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T1X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 10112001
start = start - 2.1 #  (daniel['T1_MA_02'][0]+daniel['T1_MA_03'][0])/2. - 2.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T1X']
 zpos[n+i] = zpos['T1X']+3.64+4.06+3.64-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T1X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
dt_modules['T1X'] = DtModule(tubes['T1X'],0,0,0)
for tube in dt_modules['T1X'].get_tubes():
    x = tube.get_center_position()[0]
    y = tube.get_center_position()[1]
    z = tube.get_center_position()[2]
    print("{}\t{}\t{}".format(x,y,z))
    
print("---------------------------------")

tubes['T1U'] = []
#T1u: take survey corrected points
zpos['T1U'] = (daniel['T1_MB_01'][2]+daniel['T1_MB_02'][2]+daniel['T1_MB_03'][2]+daniel['T1_MB_04'][2])/4. - 3.03 -3.64 -4.06 -3.64
angleu1 = ROOT.TMath.ATan2((daniel['T1_MB_01'][0]-daniel['T1_MB_04'][0]),(daniel['T1_MB_01'][1]-daniel['T1_MB_04'][1]))
angleu2 = ROOT.TMath.ATan2((daniel['T1_MB_02'][0]-daniel['T1_MB_03'][0]),(daniel['T1_MB_02'][1]-daniel['T1_MB_03'][1]))
angleu = (angleu1+angleu2)/2.

angle = -angleu # 60.208/180.*ROOT.TMath.Pi()  ???
#Stefan: convert angle to Euler angles
phi,theta,psi = DtAlignment.utils.z_rotation_to_euler_angles(angle)

tx,ty=0,0
for i in range(1,5):
 p = 'T1_MB_0'+str(i)
 tx += daniel[p][0] - (rn[p][0]*ROOT.TMath.Cos(angle) - rn[p][1]*ROOT.TMath.Sin(angle))
 ty += daniel[p][1] - (rn[p][0]*ROOT.TMath.Sin(angle) + rn[p][1]*ROOT.TMath.Cos(angle))

tx=tx/4.
ty=ty/4.

L =  110.
start = (rn['T1_MB_01'][0]+rn['T1_MB_04'][0])/2.
delta = (start - ( (rn['T1_MB_02'][0]+rn['T1_MB_03'][0])/2.+1.1+2.1) )/ 10.
delta = 4.2

n = 11002012
for i in range(12): 
 xnom = start - delta * i
 ynom = L/2.
 xpos[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n-i] = zpos['T1U']-deltaZ
 #build DriftTube object
 top_pos = ROOT.TVector3(xpos[n-i],ypos[n-i],zpos[n-i])
 bot_pos = ROOT.TVector3(xposb[n-i],yposb[n-i],zpos[n-i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T1U'].append(DriftTube(n-i,center[0],center[1],center[2],phi,theta,psi))
 
n = 11012001
start = (rn['T1_MB_02'][0]+rn['T1_MB_03'][0])/2.+1.1
for i in range(12): 
 xnom = start + delta * i
 ynom = L/2.
 xpos[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n+i] = zpos['T1U']+3.64-deltaZ
 top_pos = ROOT.TVector3(xpos[n+i],ypos[n+i],zpos[n+i])
 bot_pos = ROOT.TVector3(xposb[n+i],yposb[n+i],zpos[n+i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T1U'].append(DriftTube(n+i,center[0],center[1],center[2],phi,theta,psi))
n = 11102001
start = (rn['T1_MB_02'][0]+rn['T1_MB_03'][0])/2.
for i in range(12): 
 xnom = start + delta * i
 ynom = L/2.
 xpos[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n+i] = zpos['T1U']+3.64+4.06 -deltaZ
 top_pos = ROOT.TVector3(xpos[n+i],ypos[n+i],zpos[n+i])
 bot_pos = ROOT.TVector3(xposb[n+i],yposb[n+i],zpos[n+i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T1U'].append(DriftTube(n+i,center[0],center[1],center[2],phi,theta,psi))
n = 11112001
start = start - 2.1
for i in range(12): 
 xnom = start + delta * i
 ynom = L/2.
 xpos[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n+i] = zpos['T1U']+3.64+4.06+3.64-deltaZ
 top_pos = ROOT.TVector3(xpos[n+i],ypos[n+i],zpos[n+i])
 bot_pos = ROOT.TVector3(xposb[n+i],yposb[n+i],zpos[n+i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T1U'].append(DriftTube(n+i,center[0],center[1],center[2],phi,theta,psi))
 
dt_modules['T1U'] = DtModule(tubes['T1U'],0,0,0,phi,theta,psi)

# T2X:
tubes['T2X'] = []
zpos['T2X'] = (daniel['T2_MD_01'][2]+daniel['T2_MD_02'][2]+daniel['T2_MD_03'][2]+daniel['T2_MD_04'][2])/4. - 3.03 - 3.64 - 4.06 - 3.6480
ypos['T2X'] = [(daniel['T2_MD_01'][1]+daniel['T2_MD_02'][1])/2.,(daniel['T2_MD_04'][1]+daniel['T2_MD_03'][1])/2.]
n = 21112001
start = daniel['T2_MD_03'][0] # (daniel['T2_MD_02'][0]+daniel['T2_MD_03'][0])/2.  top x survey do not match nominal distance, makes 1.8mm difference to FairShip
delta = ( (daniel['T2_MD_01'][0]+daniel['T2_MD_04'][0])/2. - (start+1.1-2.1)) / 11.
delta = 4.2  # nominal distance, more reliable

for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T2X']
 zpos[n+i] = zpos['T2X']+3.64+4.06+3.64-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T2X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 21102001
start = start - 2.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T2X']
 zpos[n+i] = zpos['T2X']+3.64+4.06-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T2X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 21012001
start = start + 1.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T2X']
 zpos[n+i] = zpos['T2X']+3.64-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T2X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 21002001
start = start + 2.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 ypos[n+i] = ypos['T2X']
 zpos[n+i] = zpos['T2X']-deltaZ
 y_center = ypos[n+i][1] + ((ypos[n+i][0]-ypos[n+i][1]) / 2)
 tubes['T2X'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T2X'] = DtModule(tubes['T2X'],0,0,0)

#T2v:take survey corrected points
tubes['T2V'] = []
anglev1 = ROOT.TMath.ATan2((daniel['T2_MC_02'][0]-daniel['T2_MC_03'][0]),(daniel['T2_MC_02'][1]-daniel['T2_MC_03'][1]))
anglev2 = ROOT.TMath.ATan2((daniel['T2_MC_01'][0]-daniel['T2_MC_04'][0]),(daniel['T2_MC_01'][1]-daniel['T2_MC_04'][1]))
anglev = (anglev1+anglev2)/2.


zpos['T2V'] = (daniel['T2_MC_01'][2]+daniel['T2_MC_02'][2]+daniel['T2_MC_03'][2]+daniel['T2_MC_04'][2])/4. + 3.03
L =  110.
angle = -anglev # ???
phi, theta, psi = DtAlignment.utils.z_rotation_to_euler_angles(angle)
tx,ty=0,0
for i in range(1,5):
 p = 'T2_MC_0'+str(i)
 tx += daniel[p][0] - (rn[p][0]*ROOT.TMath.Cos(angle) - rn[p][1]*ROOT.TMath.Sin(angle))
 ty += daniel[p][1] - (rn[p][0]*ROOT.TMath.Sin(angle) + rn[p][1]*ROOT.TMath.Cos(angle))
tx=tx/4.
ty=ty/4.

n = 20112001
start = (rn['T2_MC_02'][0]+rn['T2_MC_03'][0])/2.
delta = ((rn['T2_MC_01'][0]+rn['T2_MC_04'][0])/2. - ( (rn['T2_MC_02'][0]+rn['T2_MC_03'][0])/2.-2.1+1.1)) / 11.
delta = 4.2

for i in range(12): 
 xnom = start + delta * i
 ynom = L/2.
 xpos[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n+i] = zpos['T2V']+3.64+4.06+3.64-deltaZ
 #Build tubes
 top_pos = ROOT.TVector3(xpos[n+i],ypos[n+i],zpos[n+i])
 bot_pos = ROOT.TVector3(xposb[n+i],yposb[n+i],zpos[n+i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T2V'].append(DriftTube(n+i,center[0],center[1],center[2],phi,theta,psi))
n = 20102001
start = start - 2.1
for i in range(12): 
 xnom = start + delta * i
 ynom = L/2.
 xpos[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n+i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n+i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n+i] = zpos['T2V']+3.64+4.06-deltaZ
 #Build tubes
 top_pos = ROOT.TVector3(xpos[n+i],ypos[n+i],zpos[n+i])
 bot_pos = ROOT.TVector3(xposb[n+i],yposb[n+i],zpos[n+i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T2V'].append(DriftTube(n+i,center[0],center[1],center[2],phi,theta,psi))
n = 20012012
start = (rn['T2_MC_01'][0]+rn['T2_MC_04'][0])/2.
for i in range(12): 
 xnom = start - delta * i
 ynom = L/2.
 xpos[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n-i] = zpos['T2V']+3.64-deltaZ
 #Build tubes
 top_pos = ROOT.TVector3(xpos[n-i],ypos[n-i],zpos[n-i])
 bot_pos = ROOT.TVector3(xposb[n-i],yposb[n-i],zpos[n-i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T2V'].append(DriftTube(n-i,center[0],center[1],center[2],phi,theta,psi))
n = 20002012
start = start + 2.1
for i in range(12): 
 xnom = start - delta * i
 ynom = L/2.
 xpos[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 ypos[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 ynom = -L/2.
 xposb[n-i] = xnom *ROOT.TMath.Cos(angle) - ynom*ROOT.TMath.Sin(angle) + tx
 yposb[n-i] = xnom *ROOT.TMath.Sin(angle) + ynom*ROOT.TMath.Cos(angle) + ty
 zpos[n-i] = zpos['T2V']-deltaZ
 #Build tubes
 top_pos = ROOT.TVector3(xpos[n-i],ypos[n-i],zpos[n-i])
 bot_pos = ROOT.TVector3(xposb[n-i],yposb[n-i],zpos[n-i])
 center = DtAlignment.utils.calculate_center(top_pos, bot_pos)
 tubes['T2V'].append(DriftTube(n-i,center[0],center[1],center[2],phi,theta,psi))
 
 #TODO change to real center coordinates
 dt_modules['T2V'] = DtModule(tubes['T2V'],0,0,0,phi,theta,psi)
 
#T3aX:
tubes['T3aX'] = []
zpos['T3aX'] = (( daniel['T3_T01'][2] + daniel['T3_B01'][2] + daniel['T3_T02'][2] + daniel['T3_B02'][2])/4. + 4.33)
ypos['T3aX'] = [(daniel['T3_T01'][1]+daniel['T3_T02'][1])/2.,(daniel['T3_B01'][1]+daniel['T3_B02'][1])/2.]
y_center = ypos['T3aX'][1] + ((ypos['T3aX'][0]-ypos['T3aX'][1]) / 2)

delta = ( (daniel['T3_T01'][0] + daniel['T3_B01'][0])/2. - (daniel['T3_T02'][0] + daniel['T3_B02'][0])/2. )/8.
delta = 4.2

n = 30002037
start = (daniel['T3_T02'][0] + daniel['T3_B02'][0])/2. -delta +2.1 -delta
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3aX']-deltaZ
 tubes['T3aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30012037
start = (daniel['T3_T02'][0] + daniel['T3_B02'][0])/2. -delta
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3aX']+3.64-deltaZ
 tubes['T3aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30102037
start = (daniel['T3_T02'][0] + daniel['T3_B02'][0])/2. -delta -1.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3aX']+3.64+4.06-deltaZ
 tubes['T3aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30112037
start = (daniel['T3_T02'][0] + daniel['T3_B02'][0])/2. -delta -1.1 -2.1
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3aX']+3.64+4.06+3.64-deltaZ
 tubes['T3aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T3aX'] = DtModule(tubes['T3aX'],0,0,0)

#T3bX:
tubes['T3bX'] = []
zpos['T3bX'] = (( daniel['T3_T03'][2] + daniel['T3_B03'][2] + daniel['T3_T04'][2] + daniel['T3_B04'][2])/4. + 4.33)
ypos['T3bX'] = [(daniel['T3_T03'][1]+daniel['T3_T04'][1])/2.,(daniel['T3_B03'][1]+daniel['T3_B04'][1])/2.]
y_center = ypos['T3bX'][1] + ((ypos['T3bX'][0]-ypos['T3bX'][1]) / 2)

delta = ( (daniel['T3_T03'][0] + daniel['T3_B03'][0])/2. - (daniel['T3_T04'][0] + daniel['T3_B04'][0])/2. )/8.
delta = 4.2

n = 30002025
start =  (daniel['T3_T04'][0] + daniel['T3_B04'][0])/2.  -delta +2.1 -delta
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3bX']-deltaZ
 tubes['T3bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30012025
start =  (daniel['T3_T04'][0] + daniel['T3_B04'][0])/2.  -delta 
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3bX']+3.64-deltaZ
 tubes['T3bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30102025
start =  (daniel['T3_T04'][0] + daniel['T3_B04'][0])/2.  -delta -1.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3bX']+3.64+4.06-deltaZ
 tubes['T3bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30112025
start =  (daniel['T3_T04'][0] + daniel['T3_B04'][0])/2.  -delta -1.1 -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3bX']+3.64+4.06+3.64-deltaZ
 tubes['T3bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))

#TODO change to real center coordinates
dt_modules['T3bX'] = DtModule(tubes['T3bX'],0,0,0)

#T3cX:
tubes['T3cX'] = []
zpos['T3cX'] = (( daniel['T3_T05'][2] + daniel['T3_B05'][2] + daniel['T3_T06'][2] + daniel['T3_B06'][2])/4. + 4.33)
ypos['T3cX'] = [(daniel['T3_T05'][1]+daniel['T3_T06'][1])/2.,(daniel['T3_B05'][1]+daniel['T3_B06'][1])/2.]
y_center = ypos['T3cX'][1] + ((ypos['T3cX'][0]-ypos['T3cX'][1]) / 2)

delta = ( (daniel['T3_T05'][0] + daniel['T3_B05'][0])/2. - (daniel['T3_T06'][0] + daniel['T3_B06'][0])/2. )/8.
delta = 4.2

n = 30002013
start = (daniel['T3_T06'][0] + daniel['T3_B06'][0])/2. -delta +2.1 -delta
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3cX']-deltaZ
 tubes['T3cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30012013
start = (daniel['T3_T06'][0] + daniel['T3_B06'][0])/2. -delta 
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3cX']+3.64-deltaZ
 tubes['T3cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30102013
start = start -1.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3cX']+3.64+4.06-deltaZ
 tubes['T3cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30112013
start = start -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3cX']+3.64+4.06+3.64-deltaZ
 tubes['T3cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T3cX'] = DtModule(tubes['T3cX'],0,0,0)

#T3dX:
tubes['T3dX'] = []
zpos['T3dX'] = (( daniel['T3_T07'][2] + daniel['T3_B07'][2] + daniel['T3_T08'][2] + daniel['T3_B08'][2])/4. + 4.33)
ypos['T3dX'] = [(daniel['T3_T07'][1]+daniel['T3_T08'][1])/2.,(daniel['T3_B07'][1]+daniel['T3_B08'][1])/2.]
y_center = ypos['T3dX'][1] + ((ypos['T3dX'][0]-ypos['T3dX'][1]) / 2)

delta = ( (daniel['T3_T07'][0] + daniel['T3_B07'][0])/2. - (daniel['T3_T08'][0] + daniel['T3_B08'][0])/2. )/8.
delta = 4.2


n = 30002001
start = (daniel['T3_T08'][0] + daniel['T3_B08'][0])/2. -delta +2.1 -delta
for i in range(12): 
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3dX'] -deltaZ
 tubes['T3dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30012001
start = (daniel['T3_T08'][0] + daniel['T3_B08'][0])/2. -delta 
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3dX']+3.64-deltaZ
 tubes['T3dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30102001
start =(daniel['T3_T08'][0] + daniel['T3_B08'][0])/2. -delta -1.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3dX']+3.64+4.06-deltaZ
 tubes['T3dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 30112001
start = (daniel['T3_T08'][0] + daniel['T3_B08'][0])/2. -delta -1.1 -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T3dX']+3.64+4.06+3.64-deltaZ
 tubes['T3dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T3dX'] = DtModule(tubes['T3dX'],0,0,0)

#T4aX:
tubes['T4aX'] = []
zpos['T4aX'] = (( daniel['T4_T01'][2] + daniel['T4_B01'][2] + daniel['T4_T02'][2] + daniel['T4_B02'][2])/4. - 4.33 -3.64-4.06-3.64)
ypos['T4aX'] = [(daniel['T4_T01'][1]+daniel['T4_T02'][1])/2.,(daniel['T4_B01'][1]+daniel['T4_B02'][1])/2.]
y_center = ypos['T4aX'][1] + ((ypos['T4aX'][0]-ypos['T4aX'][1]) / 2)

delta = ( (daniel['T4_T01'][0] + daniel['T4_B01'][0])/2. - (daniel['T4_T02'][0] + daniel['T4_B02'][0])/2. )/8.
delta = 4.2

n = 40002048
start = (daniel['T4_T02'][0] + daniel['T4_B02'][0])/2.  -delta +45.2
for i in range(12):
 xpos[n-i] = start - delta * i
 zpos[n-i] = zpos['T4aX']-deltaZ
 tubes['T4aX'].append(DriftTube(n-i,xpos[n-i],y_center,zpos[n-i]))
n = 40012037
start = (daniel['T4_T02'][0] + daniel['T4_B02'][0])/2.  -delta +1.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4aX']+3.64-deltaZ
 tubes['T4aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40102037
start = (daniel['T4_T02'][0] + daniel['T4_B02'][0])/2.  -delta  
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4aX']+3.64+4.06-deltaZ
 tubes['T4aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40112037
start = (daniel['T4_T02'][0] + daniel['T4_B02'][0])/2.  -delta  -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4aX']+3.64+4.06+3.64-deltaZ
 tubes['T4aX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T4aX'] = DtModule(tubes['T4aX'],0,0,0)

#T4bX:
tubes['T4bX'] = []
zpos['T4bX'] = (( daniel['T4_T03'][2] + daniel['T4_B03'][2] + daniel['T4_T04'][2] + daniel['T4_B04'][2])/4. - 4.33 -3.64-4.06-3.64)
ypos['T4bX'] = [(daniel['T4_T03'][1]+daniel['T4_T04'][1])/2.,(daniel['T4_B03'][1]+daniel['T4_B04'][1])/2.]
y_center = ypos['T4bX'][1] + ((ypos['T4bX'][0]-ypos['T4bX'][1]) / 2)

delta = ( (daniel['T4_T03'][0] + daniel['T4_B03'][0])/2. - (daniel['T4_T04'][0] + daniel['T4_B04'][0])/2. )/8.
delta = 4.2

n = 40002036
start = (daniel['T4_T04'][0] + daniel['T4_B04'][0])/2. -delta +45.2
for i in range(12):
 xpos[n-i] = start - delta * i
 zpos[n-i] = zpos['T4bX']-deltaZ
 tubes['T4bX'].append(DriftTube(n-i,xpos[n-i],y_center,zpos[n-i]))
n = 40012025
start = (daniel['T4_T04'][0] + daniel['T4_B04'][0])/2. -delta +45.2 -10*delta-2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4bX']+3.64-deltaZ
 tubes['T4bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40102025
start = (daniel['T4_T04'][0] + daniel['T4_B04'][0])/2. -delta 
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4bX']+3.64+4.06-deltaZ
 tubes['T4bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40112025
start = (daniel['T4_T04'][0] + daniel['T4_B04'][0])/2. -delta  -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4bX']+3.64+4.06+3.64-deltaZ
 tubes['T4bX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))

#TODO change to real center coordinates
dt_modules['T4bX'] = DtModule(tubes['T4bX'],0,0,0)

#T4cX:
tubes['T4cX'] = []
zpos['T4cX'] = (( daniel['T4_T05'][2] + daniel['T4_B05'][2] + daniel['T4_T06'][2] + daniel['T4_B06'][2])/4. - 4.33 -3.64-4.06-3.64)
ypos['T4cX'] = [(daniel['T4_T05'][1]+daniel['T4_T06'][1])/2.,(daniel['T4_B05'][1]+daniel['T4_B06'][1])/2.]
y_center = ypos['T4cX'][1] + ((ypos['T4cX'][0]-ypos['T4cX'][1]) / 2)

delta = ( (daniel['T4_T05'][0] + daniel['T4_B05'][0])/2. - (daniel['T4_T06'][0] + daniel['T4_B06'][0])/2. )/8.
delta = 4.2

n = 40002024
start = (daniel['T4_T06'][0] + daniel['T4_B06'][0])/2.  -delta +45.2
for i in range(12):
 xpos[n-i] = start - delta * i
 zpos[n-i] = zpos['T4cX']-deltaZ
 tubes['T4cX'].append(DriftTube(n-i,xpos[n-i],y_center,zpos[n-i]))
n = 40012013
start = (daniel['T4_T06'][0] + daniel['T4_B06'][0])/2.  -delta +45.2 -10*delta-2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4cX']+3.64-deltaZ
 tubes['T4cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40102013
start = (daniel['T4_T06'][0] + daniel['T4_B06'][0])/2.  -delta  
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4cX']+3.64+4.06-deltaZ
 tubes['T4cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40112013
start = (daniel['T4_T06'][0] + daniel['T4_B06'][0])/2.  -delta -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4cX']+3.64+4.06+3.64-deltaZ
 tubes['T4cX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T4cX'] = DtModule(tubes['T4cX'],0,0,0)

#T4dX:
tubes['T4dX'] = []
zpos['T4dX'] = (( daniel['T4_T07'][2] + daniel['T4_B07'][2] + daniel['T4_T08'][2] + daniel['T4_B08'][2])/4. - 4.33 -3.64-4.06-3.64)
ypos['T4dX'] = [(daniel['T4_T07'][1]+daniel['T4_T08'][1])/2.,(daniel['T4_B07'][1]+daniel['T4_B08'][1])/2.]
y_center = ypos['T4dX'][1] + ((ypos['T4dX'][0]-ypos['T4dX'][1]) / 2)

delta = ( (daniel['T4_T07'][0] + daniel['T4_B07'][0])/2. - (daniel['T4_T08'][0] + daniel['T4_B08'][0])/2. )/8.
delta = 4.2

n = 40002012
start = (daniel['T4_T08'][0] + daniel['T4_B08'][0])/2. -delta +45.2
for i in range(12):
 xpos[n-i] = start - delta * i
 zpos[n-i] = zpos['T4dX']-deltaZ
 tubes['T4dX'].append(DriftTube(n-i,xpos[n-i],y_center,zpos[n-i]))
n = 40012001
start = (daniel['T4_T08'][0] + daniel['T4_B08'][0])/2.  -delta +45.2 -10*delta-2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4dX']+3.64-deltaZ
 tubes['T4dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40102001
start = (daniel['T4_T08'][0] + daniel['T4_B08'][0])/2.  -delta
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4dX']+3.64+4.06-deltaZ
 tubes['T4dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
n = 40112001
start = (daniel['T4_T08'][0] + daniel['T4_B08'][0])/2.  -delta -2.1
for i in range(12):
 xpos[n+i] = start + delta * i
 zpos[n+i] = zpos['T4dX']+3.64+4.06+3.64-deltaZ
 tubes['T4dX'].append(DriftTube(n+i,xpos[n+i],y_center,zpos[n+i]))
 
#TODO change to real center coordinates
dt_modules['T4dX'] = DtModule(tubes['T4dX'],0,0,0)

for key in dt_modules.keys():
    fname = key + ".dat"
    f = open(fname,'w')
    for tube in dt_modules[key].get_tubes():
        x = tube.get_center_position()[0]
        y = tube.get_center_position()[1]
        z = tube.get_center_position()[2]
        f.write("{}\t{}\t{}\n".format(x,y,z))
    f.close()        
    
# for module in dt_modules.keys():
#     for j in range(len(dt_modules[module].get_tubes())):
#         tube = dt_modules[module].get_tubes()[j]
#         #print("Module: {}\tTube: {}\tID: {}\tLength: {}cm".format(module,j,tube._ID,tube._length))
#         print("Module: {}\tTube: {}\tLength: {}cm".format(module,j,tube._length))
    