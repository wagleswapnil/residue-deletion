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

class TransformationFE:
    def __init__(self, in_path, rerun_gmx_rerun: bool = True):
        self.in_path = in_path
        self.rerun_gmx_rerun = rerun_gmx_rerun
        self.edge1_steps = None
        self.edge3_steps = None


    def edge1_fe(self, steps):
        self.edge1_steps = steps
        self._do_edge1_sanity_check()
        edge1_fe_value, edge1_fe_error = self.calculate_edge1_fe()
        print("Edge 1 FE in kcal/mol: ")
        print(edge1_fe_value, edge1_fe_error)
        return edge1_fe_value, edge1_fe_error
    
    def edge2_fe(self):
        if self.edge1_steps == None:
            print("ERROR: No. of Edge 1 steps are needed.")
            quit()
        out_folder= os.path.join(os.getcwd(), "test", f"edge1_{self.edge1_steps}_edge3_0")
        self.calculate_reduced_potentials(f"edge1_{self.edge1_steps}", "edge3_0", out_folder)
        fe, err = calculate_fe_from_reduced_potentials(f"edge1_{self.edge1_steps}", "edge3_0", out_folder)
        print("Edge 2 FE in kcal/mol: ")
        print(fe, err)
        return fe, err
    
    def edge3_fe(self, steps):
        self.edge3_steps = steps
        subsample = False
        conservative = False
        u_nk_all_stages = []
        for stage in range(0, self.edge3_steps):
            with warnings.catch_warnings():
                # Tons of pandas warnings here
                warnings.simplefilter(action="ignore", category=FutureWarning)
                full_u_nk = extract_u_nk(
                    f"{self.in_path}/edge3_{stage}/production/production.xvg", T=298.15
                )
                u_nk = slicing(full_u_nk, step=1, lower=2000)
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
        for i in range(0, self.edge3_steps - 1):
            fe += df_value.iloc[i, i + 1]
            err += df_error.iloc[i, i + 1]

        #fe = to_kcalmol(fe, T=298.15)
        #err = to_kcalmol(err, T=298.15)

        print("Edge 3 FE in kcal/mol: ")
        print(fe, err)
        return fe, err

    def edge4_fe(self):
        if self.edge3_steps == None:
            self.edge3_steps=20
        out_folder= os.path.join(os.getcwd(), "test", f"edge3_{self.edge3_steps-1}_stage5")
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        self.calculate_reduced_potentials(f"edge3_{self.edge3_steps-1}", "stage5", out_folder)
        fe, err = calculate_fe_from_reduced_potentials(f"edge3_{self.edge3_steps-1}", "stage5", out_folder)
        print("Edge 4 FE in kcal/mol: ")
        print(fe, err)
        return fe, err


    def calculate_edge1_fe(self):
        total_fe = 0
        total_error = 0
        out_folder= os.path.join(os.getcwd(), "test", "stage1_edge1_0")
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        self.calculate_reduced_potentials("stage1", "edge1_0", out_folder)
        total_fe, total_error = calculate_fe_from_reduced_potentials("stage1", "edge1_0", out_folder)
        for step in range(0, self.edge1_steps): 
            out_folder= os.path.join(os.getcwd(), "test", f"edge1_{step}_edge1_{step + 1}")
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            self.calculate_reduced_potentials(f"edge1_{step}", f"edge1_{step + 1}", out_folder)
            fe, err = calculate_fe_from_reduced_potentials(f"edge1_{step}", f"edge1_{step + 1}", out_folder)
            total_fe += fe
            total_error += err
        return total_fe, total_error
    

    def calculate_reduced_potentials(self, state1, state2, out_folder):
        if self.rerun_gmx_rerun:
            self.run_gmx_rerun(state1, state1, out_folder)
            self._extract_potential_from_edr(state1, state1, out_folder)
        
            self.run_gmx_rerun(state1, state2, out_folder)
            self._extract_potential_from_edr(state1, state2, out_folder)
        
            self.run_gmx_rerun(state2, state1, out_folder)
            self._extract_potential_from_edr(state2, state1, out_folder)
        
            self.run_gmx_rerun(state2, state2, out_folder)
            self._extract_potential_from_edr(state2, state2, out_folder)
        return

    def _extract_potential_from_edr(self, state1, state2, out_folder):
        cmd = f"gmx energy -f {state1}_{state2}.edr -o {state1}_{state2}_energy.xvg"
        input_option = "Potential"
        subprocess.run(cmd, check=True, shell=True, text=True, input=input_option, cwd=out_folder)
        return

    def run_gmx_rerun(self, state1, state2, out_folder):
        tpr = os.path.abspath(os.path.join(self.in_path, state1, "production/production.tpr"))
        traj = os.path.abspath(os.path.join(self.in_path, state2, "production/production.xtc"))
        cmd = f"gmx mdrun -s {tpr} -rerun {traj} -deffnm {state1}_{state2}"
        subprocess.run(cmd, check=True, shell=True, cwd=out_folder)
        return

    def _do_edge1_sanity_check(self):
        for step in range(0, self.edge1_steps):
            if os.path.isdir(os.path.join(self.in_path, f"edge1_{step}")):
                if os.path.isfile(os.path.join(self.in_path, f"edge1_{step}/production", "production.xtc")) and os.path.isfile(os.path.join(self.in_path, f"edge1_{step}/production", "production.tpr")):
                    continue
                else:
                    print("production.xtc or production.tpr doesn't exist.")
                    print(os.path.join(self.in_path, f"edge1_{step}"))
                    quit()
            else:
                print("Directory doesn't exist.")
                print(os.path.join(self.in_path, f"edge1_{step}"))
                quit()
        return
