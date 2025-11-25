
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from floris import FlorisModel, TimeSeries

plt.style.use([
  "dark_background",
])

fmodel = FlorisModel("FLORIS_case.yaml")

air_density = 1.225  # kg/m^3
wind_speeds = np.arange(0, 30*1.01, 0.25)[1:]
wind_directions = 270.0*np.ones_like(wind_speeds)
wind_data = TimeSeries(
  wind_directions=wind_directions,
  wind_speeds=wind_speeds,
  turbulence_intensities=0.06,
)
wind_data.assign_ti_using_IEC_method()
fmodel.set(
  air_density=air_density,
  wind_data=wind_data,
)

diameter_rotor = 130.0  # m
area_rotor = np.pi/4*diameter_rotor**2  # m^2
fluid_force_fun = lambda V : 0.5*air_density*area_rotor*V**2
fluid_power_fun = lambda V : 0.5*air_density*area_rotor*V**3

fmodel.run()

P_T = fmodel.get_turbine_powers()/1e6  # MW
P_mean = fmodel.get_farm_power()/1e6/P_T.shape[1]  # MW

CTf_T = fmodel.get_turbine_thrust_coefficients()  # kN
Tf_T = CTf_T*np.expand_dims(fluid_force_fun(wind_speeds), -1)/1e3  # kN

fig, ax = plt.subplots()
ax.plot(wind_speeds, P_T, "-")
ax.plot(wind_speeds, P_mean, "--")

fig, ax = plt.subplots()
# ax.plot(wind_speeds, P_mean)
ax.plot(wind_speeds, Tf_T, "--")

print(wind_directions)
print(wind_speeds)
print(P_mean)

plt.show()

pd.DataFrame({
  "V": wind_speeds.tolist(),
  "P1": P_T[:,0].tolist(),
  "T1": Tf_T[:,0].tolist(),
  "P2": P_T[:,1].tolist(),
  "T2": Tf_T[:,1].tolist(),
  "P3": P_T[:,2].tolist(),
  "T3": Tf_T[:,2].tolist(),
  "Ptot": P_mean.tolist(),
  "Ttot": np.mean(Tf_T,axis=1).tolist(),
}).to_csv("floris.csv",index=False)

