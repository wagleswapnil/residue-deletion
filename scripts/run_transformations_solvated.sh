#!/bin/bash
module purge
module add gromacs/2022.1/gcc.8.4.0-cuda.11.7.1

IMPLEMENTED_REPLICATES="1 2 3"
IMPLEMENTED_PHASES="min nvt npt production"

function print_usage() {
    echo "USAGE: bash $(basename "$0") [-h] -d base_dir -x gmx -m mdp_dir -n N [-p opts] [-r replicates]"
}

if [ "$#" -eq 0 ]; then
    print_usage
    exit 1
fi

function print_help() {
    echo "Run free energy calculatoin for transformation solvated leg."
    echo
    print_usage
    echo "Options:"
    echo -e "-h\tPrint this help and exit"
    echo -e "-d dir\tThe base directory for the simulation."
    echo -e "\t The directory should contain a transfomation/solvated folder."
    echo -e "-m dir\t Directory containing all MDP files. This is the parent MDP folder."
    echo -e "-x gmx\tThe GROMACS executable to use."
    echo -e "-p \tThe run options, e.g. \"-ntmpi 2 -ntomp 8\". Optional."
    echo -e "\t No need for -maxwarn, as it will be calculated by the script."
    echo -e "-n N\tThe lambda/edgeX window to run. For a window,"
    echo -e "-r \tThe replicates to run. Default are: $IMPLEMENTED_REPLICATES"
    echo -e "\t Needed only for custom runs."
    echo
}

#  WORKDIR is the directory where the wt, mutant and transformation simulations will be run and stored.
while getopts "d:x:p:n:r:m:h" option; do
    case $option in
        d) BASEDIR="${OPTARG}" ;;
        x) GMX="${OPTARG}" ;;
        p) RUN_PARAMS="${OPTARG}" ;;
        n) RUN_NUMBER="${OPTARG}" ;;
        r) REPLICATES="${OPTARG}" ;;
        m) MDP_DIR="${OPTARG}" ;;
        h)
            print_help
            exit 0
            ;;
        *)
            print_usage
            exit 1
            ;;
    esac
done

# Make sure positional arguments are provided.
if [ -z "$BASEDIR" ] || [ -z "$GMX" ] || [ -z "$RUN_NUMBER" ] || [ -z "$MDP_DIR" ]; then
    echo "ERROR: Base directory, mdp directory, GROMACS executable, and lambda/edgeX number are required."
    print_usage
    exit 1
fi


# Check that the provided inputs are valid.
[ -d "$BASEDIR" ] || { echo "ERROR: Base directory $BASEDIR does not exist."; exit 1; }
[ -d "$MDP_DIR" ] || { echo "ERROR: MDP directory $MDP_DIR does not exist."; exit 1; }
which "$GMX" > /dev/null 2>&1 || { echo "ERROR: GROMACS executable $GMX not found in PATH."; exit 1; }

# Reassign variables.
LAMBDA_STATE=$(python3 /dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/python_workflow/run_scripts/lambda_reader.py $RUN_NUMBER)
if [ -n "$REPLICATES" ]; then
    IMPLEMENTED_REPLICATES=$REPLICATES
fi
[ -e "$BASEDIR/transformation/solvated/system.top" ] || { echo "ERROR: Topology file does not exist."; exit 1; }

# Path management.
BASEDIR=$(realpath "${BASEDIR}/transformation/solvated")
MDP_DIR=$(realpath "$MDP_DIR/transformation/solvated/fe_calcs")


# Function to execute grompp and mdrun commands.
function run_simulation(){
    local rundir=$1
    local phase=$2
    local structure=$3
    local topology=$4
    local mdp=$5
    local run_params=$6
    local replicate=$7
    local lambda_state=$8
    #NEEDED MAXWARNS: Stages1,5: 1,0. Edge1: 1,0. Edge3: 2,1
    [ "$phase" = "production" ] && warnings=2 || warnings=1

    cd "$rundir"

    if [ -f "production.xtc" ] && [ -f "production.cpt" ]; then
        MD_RUN_CMD="$GMX mdrun -deffnm $phase -cpi ${phase}.cpt ${run_params}"
        echo "Restarting simulation with command: $MD_RUN_CMD"
    else
        GROMPP_CMD="$GMX grompp -f $mdp -c $structure -p $topology -o $phase.tpr -maxwarn $warnings"
        echo "Running grompp with command: $GROMPP_CMD"
        eval $GROMPP_CMD || { echo "ERROR: grompp failed for phase $phase, replicate $replicate, lambda state $lambda_state."; exit 1; }
        MD_RUN_CMD="$GMX mdrun -deffnm $phase $run_params"
    fi
    echo "Running mdrun with command: $MD_RUN_CMD"
    eval $MD_RUN_CMD || { echo "ERROR: mdrun failed for phase $phase, replicate $replicate, lambda state $lambda_state."; exit 1; }
}

# Code to launch the simulation. 
# One simulation will be launched for a given lambda state, for both legs, all replicates, and all phases (depending on which phase needs to be run).
for REPLICATE in $IMPLEMENTED_REPLICATES; do
    echo "Running replicate $REPLICATE ..."
    WORKDIR=$BASEDIR/replicate_$REPLICATE/$LAMBDA_STATE
    TOPOLOGY=$BASEDIR/system.top
    for PHASE in $IMPLEMENTED_PHASES; do
        echo "Running $PHASE of the replicate $REPLICATE..."
        RUNDIR=$WORKDIR/$PHASE
        mkdir -p $RUNDIR
        if [ -f "$RUNDIR/${PHASE}.gro" ]; then
            continue
        fi
        
        case $PHASE in 
            "min") STRUCTURE=$BASEDIR/minimized.gro ;;
            "nvt") STRUCTURE=$WORKDIR/min/min.gro ;;
            "npt") STRUCTURE=$WORKDIR/nvt/nvt.gro ;;
            "production") STRUCTURE=$WORKDIR/npt/npt.gro ;;
            *) echo "Unknown phase $PHASE" ;;
        esac

        MDP="${MDP_DIR}/${PHASE}_${LAMBDA_STATE}.mdp"

        run_simulation "$RUNDIR" "$PHASE" "$STRUCTURE" "$TOPOLOGY" "$MDP" "$RUN_PARAMS" "$REPLICATE" "${LAMBDA_STATE}"
    done
done

            
