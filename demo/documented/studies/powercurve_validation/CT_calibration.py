import yaml

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn

################################################################################
# this is a brief script for calibration of CT.
#
# data comes from a run at 7.0 m/s (region II) with all calibration factors set
# to zero. a convergence study is run with single_turbine.yaml, and the results
# are stored and saved in CTcalibration_data.csv.
#
# we target matching the reference CT at 7.0 m/s with the asymptotic simulation.
# we assume geometric convergence of thrust w.r.t. dx. this allows the hand-fit
# of the convergence data for an estimate of the asymptotic thrust, shown below.
################################################################################

# simulation properties
rho_fluid = 1.225  # kg/m^3
diameter_rotor = 130.0  # m
area_rotor = np.pi/4*diameter_rotor**2  # m^2
V_inf = 8.0  # m/s
Tf_fluid = 0.5*rho_fluid*area_rotor*V_inf**2/1e3
P_fluid = 0.5*rho_fluid*area_rotor*V_inf**3/1e6
volume_domain = (780.0 - -390.0)*(390.0 - -390.0)*(520.0 - 0.04)

# load, amend calibration simulation data
fn_CT_data = "CTcalibration_data.csv"
df_data = pd.read_csv(fn_CT_data)
df_data["dx_elem"] = volume_domain/df_data.Nelem

# fit loge convergence with linear model
mb_Tf = np.polyfit(
  np.log10(df_data.dx_elem.to_numpy()[:-1]),
  np.log10(np.abs(rho_fluid*df_data.Tf.to_numpy()[:-1] - rho_fluid*df_data.Tf.to_numpy()[-1])),
  deg=1,
)
mb_P = np.polyfit(
  np.log10(df_data.dx_elem.to_numpy()[:-1]),
  np.log10(np.abs(rho_fluid*df_data.P.to_numpy()[:-1] - rho_fluid*df_data.P.to_numpy()[-1])),
  deg=1,
)
# i.e.:
#   log10(|J-Jinf|) = m*log(dx) + b
#   |J-Jinf| = dx**m * 10**b
#   Jinf = J -/+ dx**m * 10**b

# estimates of the true thrust based on each sample
Tf_inf_est = rho_fluid*df_data.Tf + df_data.dx_elem**mb_Tf[0] * 10**mb_Tf[1]
signpattern_P = np.ones_like(df_data.dx_elem) # np.array([1.0, 1.0, 1.0, 1.0])
P_inf_est = rho_fluid*df_data.P + signpattern_P*df_data.dx_elem**mb_P[0] * 10**mb_P[1]

# plot the convergence behavior for validation
fig, axes = plt.subplots(2, sharex=True)

# thrust vs. discretization length
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  rho_fluid*df_data.Tf.to_numpy(),
)
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  Tf_inf_est.to_numpy(),
  "--", label="corrected value",
)
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  Tf_inf_est.to_numpy()[-1]*np.ones_like(df_data.dx_elem.to_numpy()),
  "--", label="estimated asymptote",
)
axes[0].set_xlabel("discretization length, $\\Delta x$ (m)")
axes[0].set_ylabel("thrust force, $T$ (kN)")
axes[0].legend()

# apparent error vs. discretization length
axes[1].loglog(
  df_data.dx_elem.to_numpy()[:-1],
  np.abs(rho_fluid*df_data.Tf.to_numpy()[:-1] - rho_fluid*df_data.Tf.to_numpy()[-1]),
)
axes[1].loglog(
  df_data.dx_elem.to_numpy()[:-1],
  10**(mb_Tf[0]*np.log10(df_data.dx_elem.to_numpy()[:-1]) + mb_Tf[1]),
  "--", label="line of best fit",
)
axes[1].set_xlabel("discretization length, $\\Delta x$ (m)")
axes[1].set_ylabel("apparent thrust error, $e_T$ (kN)")
axes[1].legend()

# plot the convergence behavior for validation
fig, axes = plt.subplots(2, sharex=True)

# power vs. discretization length
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  rho_fluid*df_data.P.to_numpy(),
)
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  P_inf_est.to_numpy(),
  "--", label="corrected value",
)
axes[0].semilogx(
  df_data.dx_elem.to_numpy(),
  P_inf_est.to_numpy()[-1]*np.ones_like(df_data.dx_elem.to_numpy()),
  "--", label="estimated asymptote",
)
axes[0].set_xlabel("discretization length, $\\Delta x$ (m)")
axes[0].set_ylabel("power, $P$ (MW)")
axes[0].legend()

# apparent error vs. discretization length
axes[1].loglog(
  df_data.dx_elem.to_numpy()[:-1],
  np.abs(rho_fluid*df_data.P.to_numpy()[:-1] - rho_fluid*df_data.P.to_numpy()[-1]),
)
axes[1].loglog(
  df_data.dx_elem.to_numpy()[:-1],
  10**(mb_P[0]*np.log10(df_data.dx_elem.to_numpy()[:-1]) + mb_P[1]),
  "--", label="line of best fit",
)
axes[1].set_xlabel("discretization length, $\\Delta x$ (m)")
axes[1].set_ylabel("power error, $e_P$ (MW)")
axes[1].legend()

# get the reference power curve we want to match
with open("writeup/FLORIS_IEA-3p4-130-RWT.yaml", "r") as f_yaml:
  data_refPC = yaml.safe_load(f_yaml)
V_refPC = np.array(data_refPC["power_thrust_table"]["wind_speed"])
CT_refPC = np.array(data_refPC["power_thrust_table"]["thrust_coefficient"])
P_refPC = np.array(data_refPC["power_thrust_table"]["power"])/1e3

# print out the CT values we see
print(f"estimated asymptotic CT: {Tf_inf_est.to_numpy()[-1]/Tf_fluid}")
print(f"reference powercurve CT: {np.interp(V_inf, V_refPC, CT_refPC)}")
print(f"correction: {np.interp(V_inf, V_refPC, CT_refPC)/(Tf_inf_est.to_numpy()[-1]/Tf_fluid)}\n")

print(f"estimated asymptotic CP: {P_inf_est.to_numpy()[-1]/P_fluid}")
print(f"reference powercurve CP: {np.interp(V_inf, V_refPC, P_refPC)/P_fluid}")
print(f"correction: {np.interp(V_inf, V_refPC, P_refPC)/(P_inf_est.to_numpy()[-1])}")

# show and exit
plt.show()

###
