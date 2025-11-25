import pprint
import yaml

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use(
    [
        "dark_background",
        "https://raw.githubusercontent.com/cfrontin/tools_cvf/main/tools_cvf/stylesheet_cvf.mplstyle",
        "https://raw.githubusercontent.com/cfrontin/tools_cvf/main/tools_cvf/stylesheet_nrel.mplstyle",
    ]
)

with open("FLORIS_IEA-3p4-130-RWT.yaml", "r") as f_floris:
  data_FLORIS = yaml.safe_load(f_floris)

do_single_turbine=False
do_multi_turbine=True

### SINGLE TURBINE CASE

if do_single_turbine:
  df_curves_collection = []
  label_collection = []
  if True:
    # extract three levels of discretization
    casename = "kestrel_nx19ny13nz06" # "nx19ny13nz06"
    df_curves_lo = pd.read_csv(f"../curves_{casename}.csv", names=["V","P","Tf"])
    df_curves_collection.append(df_curves_lo)
    label_collection.append("$N_x=19$, $N_y=13$, $N_z=6$")
  if True:
    casename = "kestrel_nx27ny18nz09" # "nx27ny18nz09"
    df_curves_mid = pd.read_csv(f"../curves_{casename}.csv", names=["V","P","Tf"])
    df_curves_collection.append(df_curves_mid)
    label_collection.append("$N_x=27$, $N_y=18$, $N_z=9$")
  if True:
    casename = "kestrel_nx38ny25nz13" # "nx38ny25nz13"
    df_curves_hi = pd.read_csv(f"../curves_{casename}.csv", names=["V","P","Tf"])
    df_curves_collection = [df_curves_lo, df_curves_mid, df_curves_hi]
    label_collection.append("$N_x=38$, $N_y=25$, $N_z=13$")
  if True:
    casename = "kestrel_nx54ny36nz18" # "nx54ny36nz18"
    df_curves_hi = pd.read_csv(f"../curves_{casename}.csv", names=["V","P","Tf"])
    df_curves_collection = [df_curves_lo, df_curves_mid, df_curves_hi]
    label_collection.append("$N_x=54$, $N_y=36$, $N_z=18$")

  # get FLORIS reference data
  rho_fluid = 1.225
  D_rotor=130.0
  A_rotor = np.pi/4*D_rotor**2
  Ts_fluid = lambda V: 0.5*rho_fluid*A_rotor*V**2/1e3  # kN
  Ps_fluid = lambda V: 0.5*rho_fluid*A_rotor*V**3/1e6  # MW

  maxV = np.max([np.max(df.V) for df in df_curves_collection])
  V_floris = np.array(data_FLORIS["power_thrust_table"]["wind_speed"])[data_FLORIS["power_thrust_table"]["wind_speed"] <= maxV]
  P_floris = np.array(data_FLORIS["power_thrust_table"]["power"])[data_FLORIS["power_thrust_table"]["wind_speed"] <= maxV]
  Vrated_floris = 9.812675420388173
  Prated_floris = 3622.3451711828
  CPrated_floris = Prated_floris/1e3/Ps_fluid(Vrated_floris)
  CT_floris = np.array(data_FLORIS["power_thrust_table"]["thrust_coefficient"])[data_FLORIS["power_thrust_table"]["wind_speed"] <= maxV]

  ## power plot

  fig, ax = plt.subplots()

  for df_power, label_case in zip(df_curves_collection, label_collection):
    ax.plot(df_power.V, rho_fluid*df_power.P, label=label_case)
  # ax.plot(df_power_lo.V, rho_fluid*df_power_lo.P, label="$N_x=19$, $N_y=13$, $N_z=6$")
  # ax.plot(df_power_mid.V, rho_fluid*df_power_mid.P, label="$N_x=27$, $N_y=18$, $N_z=9$")
  # ax.plot(df_power_hi.V, rho_fluid*df_power_hi.P, label="$N_x=38$, $N_y=25$, $N_z=13$")
  ax.plot(
    V_floris,
    P_floris/1e3,
    "w--",
    label="IEA-3.4MW-130m reference",
  )
  ax.plot(Vrated_floris, Prated_floris/1e3, "wo")
  ax.plot(
    V_floris,
    np.minimum(
      Prated_floris/1e3*np.ones_like(V_floris),
      CPrated_floris*Ps_fluid(V_floris)),
    "w:",
    label="IEA-3.4MW-130m: $\\min(\\frac{1}{2} C_P A V^3, P_{\\mathrm{rated}}/\\rho)$"
  )
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("turbine power, $P(V)$ (MW)")
  ax.legend()
  fig.savefig("single_turb_Pconv.png", dpi=300, bbox_inches="tight")

  ## CP plot

  fig, ax = plt.subplots()

  for df_power, label_case in zip(df_curves_collection, label_collection):
    ax.plot(df_power.V, rho_fluid*df_power.P/Ps_fluid(df_power.V), label=label_case)
  # ax.plot(df_power_lo.V, rho_fluid*df_power_lo.P/Ps_fluid(df_power_lo.V), label="$N_x=19$, $N_y=13$, $N_z=6$")
  # ax.plot(df_power_mid.V, rho_fluid*df_power_mid.P/Ps_fluid(df_power_mid.V), label="$N_x=27$, $N_y=18$, $N_z=9$")
  # ax.plot(df_power_hi.V, rho_fluid*df_power_hi.P/Ps_fluid(df_power_hi.V), label="$N_x=38$, $N_y=25$, $N_z=13$")
  ax.plot(
    V_floris[1:],
    P_floris[1:]/1e3/Ps_fluid(V_floris[1:]),
    "w--",
    label="IEA-3.4MW-130m reference",
  )
  ax.plot(Vrated_floris, CPrated_floris, "wo")
  ax.plot(
    V_floris[1:],
    np.minimum(
      Prated_floris/1e3/Ps_fluid(V_floris[1:]),
      CPrated_floris),
    "w:",
    label="IEA-3.4MW-130m: $\\min(C_P, P_{\\mathrm{rated}}/P_{\\mathrm{fluid}})$"
  )
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("power_coefficient, $C_P(V)$ (-)")
  ax.legend()
  fig.savefig("single_turb_CPconv.png", dpi=300, bbox_inches="tight")

  # fig, ax = plt.subplots()
  #
  # for df_power, label_case in zip(df_curves_collection, label_collection):
  #   ax.plot(
  #     df_power.V,
  #     (
  #       rho_fluid*df_power.P/Ps_fluid(df_power.V)
  #     )/np.interp(
  #       df_power.V,
  #       V_floris[1:],
  #       P_floris[1:]/1e3/Ps_fluid(V_floris[1:]),
  #     ),
  #     label=label_case,
  #   )

  ## thrust plot

  fig, ax = plt.subplots()

  for df_thrust, label_case in zip(df_curves_collection, label_collection):
    ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf, label=label_case)
  # ax.plot(df_thrust_lo.V, rho_fluid*df_thrust_lo.Tf, label="$N_x=19$, $N_y=13$, $N_z=6$")
  # ax.plot(df_thrust_mid.V, rho_fluid*df_thrust_mid.Tf, label="$N_x=27$, $N_y=18$, $N_z=9$")
  # ax.plot(df_thrust_hi.V, rho_fluid*df_thrust_hi.Tf, label="$N_x=38$, $N_y=25$, $N_z=13$")
  ax.plot(
    V_floris,
    CT_floris*Ts_fluid(V_floris),
    "w--",
    label="IEA-3.4MW-130m reference",
  )
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("turbine thrust, $T(V)$ (kN)")
  ax.legend()
  fig.savefig("single_turb_Tconv.png", dpi=300, bbox_inches="tight")

  ## CT plot

  fig, ax = plt.subplots()

  for df_thrust, label_case in zip(df_curves_collection, label_collection):
    ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf/Ts_fluid(df_thrust.V), label=label_case)
  # ax.plot(df_thrust_lo.V, rho_fluid*df_thrust_lo.Tf/Ts_fluid(df_thrust_lo.V), label="$N_x=19$, $N_y=13$, $N_z=6$")
  # ax.plot(df_thrust_mid.V, rho_fluid*df_thrust_mid.Tf/Ts_fluid(df_thrust_mid.V), label="$N_x=27$, $N_y=18$, $N_z=9$")
  # ax.plot(df_thrust_hi.V, rho_fluid*df_thrust_hi.Tf/Ts_fluid(df_thrust_hi.V), label="$N_x=38$, $N_y=25$, $N_z=13$")
  ax.plot(
    V_floris,
    CT_floris,
    "w--",
    label="IEA-3.4MW-130m reference",
  )
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("thrust coefficient, $C_T(V)$ (-)")
  ax.legend()
  fig.savefig("single_turb_CTconv.png", dpi=300, bbox_inches="tight")



### MULTI TURBINE CASE

if do_multi_turbine:
  df_curves_collection = []
  label_collection = []
  if True:
    # extract three levels of discretization
    casename = "kestrel_nx64ny13nz06" # "nx64ny13nz06"
    df_curves_lo = pd.read_csv(f"../curves_multi_{casename}.csv", names=["V","P1","Tf1","P2","Tf2","P3","Tf3","P_tot","Tf_tot"])
    df_curves_collection.append(df_curves_lo)
    label_collection.append("$N_x=64$, $N_y=13$, $N_z=6$")
  if True:
    # extract three levels of discretization
    casename = "kestrel_nx90ny18nz09" # "nx90ny18nz09"
    df_curves_lo = pd.read_csv(f"../curves_multi_{casename}.csv", names=["V","P1","Tf1","P2","Tf2","P3","Tf3","P_tot","Tf_tot"])
    df_curves_collection.append(df_curves_lo)
    label_collection.append("$N_x=90$, $N_y=18$, $N_z=9$")
  if True:
    # extract three levels of discretization
    casename = "kestrel_nx127ny25nz13" # "nx127ny25nz13"
    df_curves_lo = pd.read_csv(f"../curves_multi_{casename}.csv", names=["V","P1","Tf1","P2","Tf2","P3","Tf3","P_tot","Tf_tot"])
    df_curves_collection.append(df_curves_lo)
    label_collection.append("$N_x=127$, $N_y=25$, $N_z=13$")
  if False:
    # extract three levels of discretization
    casename = "kestrel_nx180ny36nz18" # "nx180ny36nz18"
    df_curves_lo = pd.read_csv(f"../curves_multi_{casename}.csv", names=["V","P1","Tf1","P2","Tf2","P3","Tf3","P_tot","Tf_tot"])
    df_curves_collection.append(df_curves_lo)
    label_collection.append("$N_x=180$, $N_y=36$, $N_z=18$")

  # get FLORIS reference data
  rho_fluid = 1.225
  D_rotor=130.0
  A_rotor = np.pi/4*D_rotor**2
  Ts_fluid = lambda V: 0.5*rho_fluid*A_rotor*V**2/1e3  # kN
  Ps_fluid = lambda V: 0.5*rho_fluid*A_rotor*V**3/1e6  # MW

  df_floris_multi = pd.read_csv("floris.csv")
  V_floris = df_floris_multi.V
  P1_floris = df_floris_multi.P1
  P2_floris = df_floris_multi.P2
  P3_floris = df_floris_multi.P3
  Tf1_floris = df_floris_multi.T1
  Tf2_floris = df_floris_multi.T2
  Tf3_floris = df_floris_multi.T3

  ## power plot

  fig, ax = plt.subplots()

  ax.plot([], [], "w--", label="T1")
  ax.plot([], [], "w:", label="T2")
  ax.plot([], [], "w-.", label="T3")

  for df_power, label_case in zip(df_curves_collection, label_collection):
    pt0 = ax.plot([], [], "-", label=label_case)
    ax.plot(df_power.V, rho_fluid*df_power.P1, "--", c=pt0[-1].get_color(), label="__")
    ax.plot(df_power.V, rho_fluid*df_power.P2, ":", c=pt0[-1].get_color(), label="__")
    ax.plot(df_power.V, rho_fluid*df_power.P3, "-.", c=pt0[-1].get_color(), label="__")
    # ax.plot(df_power.V, rho_fluid*df_power.P_tot/3, "-", c=pt0[-1].get_color(), label=label_case)
  ax.plot([], [], "w-", label="FLORIS")
  ax.plot(V_floris, P1_floris, "w--", label="_FLORIS T1_")
  ax.plot(V_floris, P2_floris, "w:", label="_FLORIS T2_")
  ax.plot(V_floris, P3_floris, "w-.", label="_FLORIS T3_")
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("turbine power, $P(V)$ (MW)")
  ax.legend()
  fig.savefig("multi_turb_Pconv.png", dpi=300, bbox_inches="tight")

  ## thrust plot

  fig, ax = plt.subplots()

  ax.plot([], [], "w--", label="T1")
  ax.plot([], [], "w:", label="T2")
  ax.plot([], [], "w-.", label="T3")

  for df_thrust, label_case in zip(df_curves_collection, label_collection):
    pt0 = ax.plot([], [], "-", label=label_case)
    ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf1, "--", c=pt0[-1].get_color(), label="__")
    ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf2, ":", c=pt0[-1].get_color(), label="__")
    ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf3, "-.", c=pt0[-1].get_color(), label="__")
    # ax.plot(df_thrust.V, rho_fluid*df_thrust.Tf_tot/3, "-", c=pt0[-1].get_color(), label=label_case)
  ax.plot([], [], "w-", label="FLORIS")
  ax.plot(V_floris, Tf1_floris, "w--", label="_FLORIS T1_")
  ax.plot(V_floris, Tf2_floris, "w:", label="_FLORIS T2_")
  ax.plot(V_floris, Tf3_floris, "w-.", label="_FLORIS T3_")
  ax.set_xlabel("hub-height farfield velocity, $V$ (m/s)")
  ax.set_ylabel("turbine thrust, $T(V)$ (kN)")
  ax.legend()
  fig.savefig("multi_turb_Tconv.png", dpi=300, bbox_inches="tight")



plt.show()
