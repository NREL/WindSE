
# WindSE Power Curve Validation

`WindSE` has the ability to feather thrust and power coefficients in order to emulate a controlled turbine through multiple regions of the power curve.
In this study, we will validate the power curve behavior in isolation and with multiple inline turbines.
We consider the discretization dependence of the `WindSE` results, and the comparison of `WindSE` power estimates compared to other tools.

## Tuning methodology

Using the `IEA-3.4MW-130m` turbine, we run a calibration study at one region II velocity and one region III velocity to adjust the power curve model parameters to match the reference curves for the turbine.
A magic number scaling parameter is applied to the output power to achieve a match between the force/velocity product power and the expected output shaft power.

## Single turbine: discretization dependence

In this section we compare the `WindSE` results to the canonical IEA 3.4MW 130m turbine's power curve that is an input to `FLORIS`.

![Power curve convergence plot](single_turb_Tconv.png)

The first plot here shows the thrust curves compared to the reference turbine data, and demonstrates a strong match for a single turbine.
The thrust peaking behavior in the reference power curve is due to the lack of thrust peak-shaving control, and `WindSE` implicitly performs peak shaving by the smoothing method used to achieve differentiability.
Long run behavior in velocity exhibits some inexact scaling; this is likely due to inaccuracy-- or at least mismatch-- in the induction profiles assumed for `WindSE`'s actuator disk.

![Power curve convergence plot](single_turb_Pconv.png)

FLORIS includes cut-in behavior that is not modeled in region I or region I.5.
Comparison with the standard IEA-3.4MW-130m power curve used for FLORIS demonstrates similar power behavior in region II, with an slight apparent overprediction of $C_P$.
The FLORIS curve does not have thrust-peak shaving, while the `WindSE` power curve implicitly does (via the smoothing), so there is a power mismatch in region II.5.

![$C_P$ curve convergence plot](single_turb_CPconv.png)

Plotting the power coefficient indicates that $C_P$ is in fact initially overpredicted in region II, and moreover $C_P$ is not quite constant in region II.
As noted previously, region II.5 is not modeled in the FLORIS power curve, but long run behavior of $C_P$ appears to be accurate, in region III.

## Multi-turbine: discretization dependence

We now run a three-turbine case with waked turbines to validate that behavior matches expectation.
In this case, we have three directly-waked turbines at $7.5D$ spacing.

![Farm thrust curve convergence plot](multi_turb_Tconv.png)

As previously, FLORIS does not implement peak-shaving which is implicitly included in the WindSE results.
Outside of region II.5 (peak-shaving), there is a strong qualitative match between the thrust of the turbines in WindSE and the FLORIS simulation.
This indicates that WindSE is correctly estimating the thrust on a given turbine as well as passing the resulting thrust into the flowfield correctly.

![Farm power curve convergence plot](multi_turb_Pconv.png)

In the power curves, WindSE again makes a strong prediction that is closely matched with the FLORIS result, less thrust peak-shaving in region II.5.


<!-- -->
