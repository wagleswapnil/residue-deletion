import subprocess
import os
import alchemlyb
from resdel.analysis.calculate_fe_from_reduced_potentials import calculate_fe_from_reduced_potentials
from glob import glob
import warnings
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.postprocessors.units import get_unit_converter
from alchemlyb.preprocessing import slicing, statistical_inefficiency
from alchemlyb.preprocessing.subsampling import decorrelate_u_nk
from alchemlyb.visualisation import plot_mbar_overlap_matrix
from alchemlyb.postprocessors.units import to_kcalmol

class AbsoluteFE:
    def __init__(self, in_path):
        self.in_path = in_path
        self.lambda_steps = None

    def total_fe(self, steps):
        self.lambda_steps = steps
        subsample = False
        conservative = True
        u_nk_all_stages = []
        for stage in range(0, self.lambda_steps):
            with warnings.catch_warnings():
                # Tons of pandas warnings here
                warnings.simplefilter(action="ignore", category=FutureWarning)
                full_u_nk = extract_u_nk(
                    f"{self.in_path}/lambda_{stage}/production/production.xvg", T=298.15
                )
                u_nk = slicing(full_u_nk, step=1)
                if subsample:
                    u_nk_all_stages.append(
                        statistical_inefficiency(
                            u_nk,
                            series=u_nk[u_nk.columns[stage]],
                            conservative=conservative,
                        )
                    )
                else:
                    u_nk_all_stages.append(u_nk)
        decorrelated_u_nk_list = [decorrelate_u_nk(u_nk) for u_nk in u_nk_all_stages]
        mbar=MBAR().fit(alchemlyb.concat(decorrelated_u_nk_list))
        #plot_mbar_overlap_matrix(mbar.overlap_matrix)
        df_value = to_kcalmol(mbar.delta_f_, T=298.15)
        df_error = to_kcalmol(mbar.d_delta_f_, T=298.15)
        fe = 0.0
        err = 0.0
        for i in range(0, self.lambda_steps - 1):
            fe += df_value.iloc[i, i + 1]
            err += df_error.iloc[i, i + 1]

        #fe = to_kcalmol(fe, T=298.15)
        #err = to_kcalmol(err, T=298.15)

        print("Total FE in kcal/mol: ")
        print(fe, err)
        return fe, err
