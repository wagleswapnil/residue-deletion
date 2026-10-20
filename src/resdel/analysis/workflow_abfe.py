import sys
from resdel.analysis.abfe_fe_calculator import AbsoluteFE

in_path = sys.argv[1]
free_energy_calculator = AbsoluteFE(in_path)
total_fe_value, total_fe_error = free_energy_calculator.total_fe(steps=20) #ToDo-- Should it be 15 or 16? Make it consistent with Edge 3 call 2 lines below.




#tpr1 = sys.argv[1]
#traj1 = sys.argv[2]
#tpr2 = sys.argv[3]
#traj2 = sys.argv[4]

#run_gmx_rerun(tpr1, traj1, tpr2, traj2)
#calculate_reduced_potential(tpr, traj)
