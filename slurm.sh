#! /bin/bash

# Submits job to a Slurm job scheduler.
# (see https://slurm.schedmd.com/sbatch.html)

# SCRIPT DIRECTORY

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"  # directory where this script is

# DEFAULT PARAMETERS
# We are concerned here with the `xmaris' computer of the Instituut-Lorentz
# at Universiteit Leiden.
# (see https://helpdesk.lorentz.leidenuniv.nl/wiki/doku.php?id=institute_lorentz:xmaris)

_EXE_DIR=/marisdata/$LOGNAME    # default script execution directory
_OUT_DIR=${_EXE_DIR}/out        # default standard output directory
_ERR_DIR=${_EXE_DIR}/err        # default standard error output directory

_PARTITION=compAMD  # default partition for the ressource allocation
_GRES=              # default generic consumable ressources
_NODES=1            # default required number of nodes
_NTASKS=1           # default number of MPI ranks running per node
_ARRAY_SIZE=        # default job array size
_ARRAY_TASKS=       # default maximum number of simultaneous tasks in the array
_TIME=              # default required time
_MEMORY=            # default real memory required per node

# HELP MENU

usage() {

less <<< "Submit job to a Slurm job scheduler.
(see https://slurm.schedmd.com/sbatch.html)
(see https://helpdesk.lorentz.leidenuniv.nl/wiki/doku.php?id=institute_lorentz:xmaris)

SYNOPSIS

  [bash] slurm.sh [OPTIONS] [ENVIRONMENT VARIABLES] [SCRIPT]

OPTIONS

  -h    Display this help.

  -w    Pause this script until completion of the job.

  -j    Job name on Slurm scheduler.
        DEFAULT: script name after last '/'
  -c    Execute after job with this ID has succesfully executed.
        DEFAULT: (not specified)

  -d    Directory in which to execute the script.
        DEFAULT: $_EXE_DIR
  -o    Standard output directory.
        NOTE: Output files are named according to job ID.
        DEFAULT: $_OUT_DIR
  -e    Standard error output directory.
        NOTE: Error files are named according to job ID.
        DEFAULT: $_ERR_DIR

  -p    Partition for the resource allocation.
        DEFAULT: $_PARTITION
  -g    Generic consumable resources.
        DEFAULT: $_GRES
  -n    Required number of nodes.
        DEFAULT: $_NODES
  -r    Number of MPI ranks running per node.
        Number of threads for OpenMP parallelised jobs.
        DEFAULT: $_NTASKS
  -a    Job array size.
        (see https://slurm.schedmd.com/job_array.html)
        NOTE: SLURM_ARRAY_TASK_ID is set as task id (between 0 and size - 1).
        DEFAULT: $_ARRAY_SIZE
  -s    Maximum number of simultaneously running tasks in the job array.
        (see https://slurm.schedmd.com/job_array.html)
        NOTE: An empty string will not set this maximum.
        DEFAULT: $_ARRAY_TASKS
  -t    Required time.
        DEFAULT: $_TIME
  -m    Real memory required per node.
        NOTE: MaxMemPerNode allocates maximum memory.
        DEFAULT: $_MEMORY
"
}

# OPTIONS

while getopts "hwj:c:d:o:e:p:g:n:r:a:s:t:m:" OPTION; do
    case $OPTION in

    h)  # help menu
        usage; exit 0;;

    w)  # wait
        WAIT=true;;

    j)  # job name
        JOB_NAME=$OPTARG;;
    c)  # chained job
        CHAIN=$OPTARG;;

    d)  # script execution directory
        EXE_DIR=$OPTARG;;
    o)  # standard output directory
        OUT_DIR=$OPTARG;;
    e)  # standard error output directory
        ERR_DIR=$OPTARG;;

    p)  # partition
        PARTITION=$OPTARG;;
    g)  # generic consumable resources
        GRES=$OPTARG;;
    n)  # nodes
        NODES=$OPTARG;;
    r)  # tasks
        NTASKS=$OPTARG;;
    a)  # array size
        ARRAY_SIZE=$OPTARG;;
    s)  # array tasks
        ARRAY_TASKS=$OPTARG;;
    t)  # time
        TIME=$OPTARG;;
    m)  # real memory
    MEMORY=$OPTARG;;

esac
done
shift $(expr $OPTIND - 1);

if [[ -z "$@" ]]; then
  echo 'No script submitted.';
  usage;
  exit 1;
fi

SCRIPT=$@   # script to execute

# JOB PARAMETERS

JOB_NAME=${JOB_NAME-${SCRIPT##*/}}  # job name

EXE_DIR=${EXE_DIR-$_EXE_DIR}; mkdir -p "$EXE_DIR";  # script execution directory
OUT_DIR=${OUT_DIR-$_OUT_DIR}; mkdir -p "$OUT_DIR";  # standard output directory
ERR_DIR=${ERR_DIR-$_ERR_DIR}; mkdir -p "$ERR_DIR";  # standard error output directory

PARTITION=${PARTITION-$_PARTITION}        # partition for the resource allocation
GRES=${GRES-$_GRES}                       # generic consumable resources
NODES=${NODES-$_NODES}                    # required number of nodes
NTASKS=${NTASKS-$_NTASKS}                 # maximum ntasks to be invoked on each core
ARRAY_SIZE=${ARRAY_SIZE-$_ARRAY_SIZE}     # job array size
ARRAY_TASKS=${ARRAY_TASKS-$_ARRAY_TASKS}  # maximum number of simultaneous tasks in the array
TIME=${TIME-$_TIME}                       # required time
MEMORY=${MEMORY-$_MEMORY}                 # real memory required per node

# SUBMIT JOB

sbatch ${WAIT:+-W} ${CHAIN:+-d afterok:$CHAIN} <<EOF
#! /bin/bash
#SBATCH --job-name='$JOB_NAME'
#SBATCH -D '$EXE_DIR'
#SBATCH --output='${OUT_DIR}/%j'
#SBATCH --error='${ERR_DIR}/%j'
#SBATCH --partition=$PARTITION
${GRES:+#SBATCH --gres=$GRES}
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NTASKS
${ARRAY_SIZE:+#SBATCH --array=0-$(($ARRAY_SIZE-1))${ARRAY_TASKS+%$ARRAY_TASKS}}
${TIME:+#SBATCH --time=$TIME}
${MEMORY:+#SBATCH --mem=$MEMORY}

export OMP_NUM_THREADS=$NTASKS

# print job parameters to both standard output and standard error output files
(>&1 printf '%-21s: %s\n' 'SUBMIT DIRECTORY' '$(pwd)')
(>&1 printf '%-21s: %s\n' 'DATE' '$(date)')
(>&1 echo)
(>&1 printf '%-21s: %s\n' 'JOB NAME' '$JOB_NAME')
(>&1 echo)
(>&1 printf '%-21s: %s\n' 'SIMULATION DIRECTORY' '$SIM_DIR')
(>&1 printf '%-21s: %s\n' 'OUTPUT FILE' '$OUT_FILE')
(>&1 echo)
(>&1 printf '%-21s: %s\n' 'PARTITION' '$PARTITION')
(>&1 printf '%-21s: %s\n' 'GRES' '$GRES')
(>&1 printf '%-21s: %s\n' 'NODES REQUIRED' '$NODES')
(>&1 printf '%-21s: %s\n' 'TASKS PER NODE' '$NTASKS')
(>&1 printf '%-21s: %s\n' 'ARRAY SIZE' '$ARRAY_SIZE')
(>&1 printf '%-21s: %s\n' 'TASKS IN ARRAY' '$ARRAY_TASKS')
(>&1 printf '%-21s: %s\n' 'TIME REQUIRED' '$TIME')
(>&1 printf '%-21s: %s\n' 'MEMORY REQUIRED' '$MEMORY')
(>&1 echo)
(>&1 printf '%-21s: %s\n' 'SCRIPT' '$SCRIPT')
(>&1 echo)

module load Python/3.10.4-GCCcore-11.3.0                    # load python3.10

$SCRIPT                                                     # launching script
EOF

${WAIT:+wait}   # wait until completion of the job

