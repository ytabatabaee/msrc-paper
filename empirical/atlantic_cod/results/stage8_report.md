# Atlantic cod Stage-8 report

Stage 8 is a post-freeze projection/sensitivity analysis. It uses Stage-7 SNAPP MRCA times and fixed Stage-6 ASTRAL topologies, and it does not rerun SNAPP or ASTRAL4.

## Primary T0_all summary

Selected root edge: `Gadmor_bat_spc+Gadmor_twc_spc+Gadmor_two_spc|Gadmor_avc_spc+Gadmor_avo_spc+Gadmor_bor_spc+Gadmor_icc_spc+Gadmor_ico_spc+Gadmor_kie_spc+Gadmor_lfc_spc+Gadmor_lfo_spc+Gadmor_low_spc`.
Root age outside/all/shift: 0.0590363 / 0.0656562 / 0.00661989.
Median absolute internal-node shift: 0.00413352.
Maximum absolute node-age shift: 0.012448.
Mean signed node-age shift: 0.0010419.
Older/younger/unchanged nodes with all windows: 6 / 5 / 0.

## Projection fit

| topology | treatment | objective | RMSE | MAE | correlation | max abs residual |
|---|---|---|---:|---:|---:|---:|
| T0_all | outside | L2 | 0.0012225 | 0.000889311 | 0.995607 | 0.00387704 |
| T0_all | all | L2 | 0.00652975 | 0.00486227 | 0.913166 | 0.0150307 |
| T0_all | outside_L1_approx | L1_median_target_approx | 0.00124733 | 0.000897821 | 0.995505 | 0.00348163 |
| T0_all | all_L1_approx | L1_median_target_approx | 0.00663442 | 0.00478176 | 0.911881 | 0.0167296 |
| T1_outside | outside | L2 | 0.00106975 | 0.000811561 | 0.996638 | 0.00287735 |
| T1_outside | all | L2 | 0.00759219 | 0.00618814 | 0.880575 | 0.0168831 |
| T1_outside | outside_L1_approx | L1_median_target_approx | 0.00112576 | 0.000785951 | 0.996446 | 0.00353081 |
| T1_outside | all_L1_approx | L1_median_target_approx | 0.00789428 | 0.00599336 | 0.872557 | 0.0204495 |

## Secondary T1_outside summary

Root age outside/all/shift: 0.0590363 / 0.0656562 / 0.00661989.
Median absolute internal-node shift: 0.00241876.
