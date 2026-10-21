import sys
from resdel.analysis.transformation_fe_calculator import TransformationFE


def run_analyze_workflow(in_path):
    #in_path = sys.argv[1]
    free_energy_calculator = TransformationFE(in_path, rerun_gmx_rerun=False)
    edge1_fe_value, edge1_fe_error = free_energy_calculator.edge1_fe(steps=8) #ToDo-- Should it be 15 or 16? Make it consistent with Edge 3 call 2 lines below.
    edge2_fe_value, edge2_fe_error = free_energy_calculator.edge2_fe()
    #edge3_fe_value, edge3_fe_error = free_energy_calculator.edge3_fe(steps=20)
    edge4_fe_value, edge4_fe_error = free_energy_calculator.edge4_fe()
    return
