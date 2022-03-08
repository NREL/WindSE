import numpy as np
import os

# Path to openfast AD15 file of the OpenFAST model to be imported into WindSE
path2ad15 = '/Users/pbortolo/work/2_openfast/r-test/glue-codes/openfast/5MW_Land_BD_DLL_WTurb/NRELOffshrBsline5MW_Onshore_AeroDyn15.dat'
destination_folder = '/Users/pbortolo/work/5_windse/WindSE/demo/documented/Yaml_Examples/Input_Data/nrel_5mw'

if os.path.isfile(path2ad15):
    f = open(path2ad15)
else:
    raise Exception('The path specified does not point to any file. Please check the path specified in line 5')

of_folder = os.path.dirname(path2ad15)

# Build the paths to airfoil polar files and blade aero data
while 1:
    string = f.readline()
    if string.split()[1] == 'NumAFfiles':
        NumAFfiles = int(string.split()[0])
        path2af = [''] * NumAFfiles
        path2af[0] = os.path.join(of_folder, f.readline().split()[0].replace('"',''))
        for i in range(1,NumAFfiles):
            path2af[i] = os.path.join(of_folder, f.readline()[:-2].replace('"',''))
        f.readline()
        f.readline()
        ad15blfilename = f.readline().split()[0].replace('"','')
        path2ad15bl = os.path.join(of_folder, ad15blfilename)
        f.close()
        break

# Load and write blade data
ad15_blade_data = np.loadtxt(path2ad15bl, skiprows=6)
windse_blade_data = np.zeros((len(ad15_blade_data[:,0]), 9))
# Non dimensional blade span (-) (ad15 has dimensional blade span (m))
windse_blade_data[:,0] = ad15_blade_data[:,0] / ad15_blade_data[-1,0]
# Chord (m)
windse_blade_data[:,1] = ad15_blade_data[:,5]
# Twist (deg)
windse_blade_data[:,2] = ad15_blade_data[:,4]
np.savetxt(os.path.join(destination_folder, 'blade_data.csv'), windse_blade_data, delimiter=",",
    header='Blade n.d. span [-], 	 Chord [m], 	 Twist [deg], 	 Lift coeff [-], 	 Drag coeff [-], 	 Induction [-], 	 Angle of attack [deg], 	 Lift force [N/m], 	 Drag force [N/m]')

# Load and write airfoil polars
af_index = np.array(ad15_blade_data[:,-1], dtype=int)
windse_path2af = os.path.join(destination_folder, 'airfoil_polars')
if not os.path.isdir(windse_path2af):
    os.makedirs(windse_path2af)
for i in range(len(af_index)):
    f = open(path2af[af_index[i]-1])
    while 1:
        string = f.readline()
        if len(string.split()) > 1:
            if string.split()[1] == 'NumTabs':
                if int(string.split()[0]) > 1:
                    raise Warning('Only the first polar table is read, but multiple are provided in the airfoil file ' + path2af[i] +  ' from openFAST!')
            if string.split()[1] == 'NumAlf':
                NumAlf = int(string.split()[0])
                polars = np.zeros((NumAlf, 4))
                f.readline()
                f.readline()
                for j in range(NumAlf):
                    polars[j,:] = np.array(f.readline().split(), dtype=float)
                polars[:,0] = np.deg2rad(polars[:,0])
                f.close()
                break
    np.savetxt(os.path.join(windse_path2af, 
        'af_station_' + str(i) + '.txt'),polars, 
        header='Angle of attack [rad] 	 Lift coefficient [-] 	 Drag coefficient [-] 	 Moment coefficient [-]')