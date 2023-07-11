import numpy as np
import matplotlib.pyplot as plt
from weis.aeroelasticse.FAST_post import FAST_IO_timeseries

outb_filename = '10ms_nss.outb'
plot = True
destination_csv = 'WindSE/demo/documented/Yaml_Examples/Input_Data/iea_15/ptfm_motion_10ms_nss.csv'

outb_data = FAST_IO_timeseries(outb_filename)
csv_data = np.zeros((len(outb_data['Time']), 7))
csv_data[:,0] = outb_data['Time'] - outb_data['Time'][0] # Start at 0 seconds
csv_data[:,1] = outb_data['PtfmSurge']
csv_data[:,2] = outb_data['PtfmSway']
csv_data[:,3] = outb_data['PtfmHeave']
csv_data[:,4] = outb_data['PtfmRoll']
csv_data[:,5] = outb_data['PtfmPitch']
csv_data[:,6] = outb_data['PtfmYaw']

np.savetxt(destination_csv, csv_data, header = 'Time (s) \t PtfmSurge (m) \t PtfmSway (m) \t PtfmHeave (m) \t PtfmRoll (deg) \t PtfmPitch (deg) \t PtfmYaw (deg)')

if plot:
    f, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].plot(csv_data[:,0], csv_data[:,1], label = 'Surge (m)')
    ax[0].plot(csv_data[:,0], csv_data[:,2], label = 'Sway (m)')
    ax[0].plot(csv_data[:,0], csv_data[:,3], label = 'Heave (m)')
    ax[1].plot(csv_data[:,0], csv_data[:,4], label = 'Roll (deg)')
    ax[1].plot(csv_data[:,0], csv_data[:,5], label = 'Pitch (deg)')
    ax[1].plot(csv_data[:,0], csv_data[:,6], label = 'Yaw (deg)')
    ax[0].set_xlabel('Time (s)', fontweight = 'bold')
    ax[0].set_ylabel('Displacement (m)', fontweight = 'bold')
    ax[1].set_xlabel('Time (s)', fontweight = 'bold')
    ax[1].set_ylabel('Rotation (deg)', fontweight = 'bold')
    ax[0].legend()
    ax[1].legend()
    ax[0].grid(color=[0.8, 0.8, 0.8], linestyle="--")
    ax[1].grid(color=[0.8, 0.8, 0.8], linestyle="--")
    plt.tight_layout()
    plt.savefig('10ms_nss')
    plt.show()