#!/bin/bash
set -x
set -euo pipefail
# Parcourt tous les dossiers h_*
for dir in 2025_07_24/h_*; do
    [ -d "$dir" ] || continue  # skip non-dossiers
    name=$(basename "$dir")
    raw="${name#h_}"

    float_value=$(julia --quiet --eval "
    try
        x = parse(Float64, \"$raw\")
        println(x)
    catch
        exit(1)
    end
    ")
    # Vérifie que la variable a été remplie
    if [[ -z "$float_value" ]]; then
       echo "$raw n'est pas un Float64 valide pour Julia"
       exit 1
    fi
    echo "$name  Float64 validé = $float_value"

    logfile="$log_${name}.log"                  # ex: log_h_0.513.log
    output_json="${dir}/output_${name}.json"          # ex: output_h_0.513.json
    input_template="${dir}/input_${name}.json"        # fichier avec placeholder
    input_real="${dir}/tmp_input_${name}.json"        # fichier reel avec job_id
    slurm_file="${dir}/job_${name}.slurm"      # fichier .slurm genereé

    # ----- 1. Cree le fichier JSON modele -----
    cat > "${input_template}" <<EOF
{
 	"job_id" : "TO_BE_REPLACED_BY_SLURM",
	"parameter" : "${name}",
	"step beta" : 0.1,
        "N" : 100,
        "J" : 1,
        "maximum value for beta" : 10,
        "step value for beta list" : 0.1,
        "cutoff" : 1e-15,
        "max bond dimension" :  300,
        "Trotter-Suzuki step" : 1e-3,
        "fixed disorder" : $float_value, 
        "gammescale" : 0.8,
        "axis" : "z",
        "savefile" : "${output_json}"
}
EQF
     jq empty "${input_template}" || { echo "JSON invalide dans $input_template"; exit 1; }

     abs_dir=$(realpath "$dir")

     echo $abs_dir

    # ----- 2. Creer le fichier SLURM -----
    cat > "${slurm_file}" <<EOF
#!/bin/bash
#SBATCH --partition=genx
#SBATCH --constraint=ib-icelake
#SBATCH --cpus-per-task=1
#SBATCH --mem 10G
#SBATCH --job-name=$name
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --threads-per-core=1    # 2 => hyperthreading
#SBATCH --time 0-72:0:0
#SBATCH --output=$logfile
#SBATCH --error=$logfile
#SBATCH --chdir=$abs_dir

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=false

echo STARTING AT $(date '+%A %d %B %Y %H:%M')
echo HOST $(hostname)
echo JOB ID $SLURM_JOB_ID
echo PARTITION $SLURM_JOB_PARTITION
echo QOS $SLURM_JOB_QOS
echo MEMORY $SLURM_MEM_PER_NODE
echo SLURM_CPUS_ON_NODE $SLURM_CPUS_ON_NODE
echo SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK
echo batch file $(squeue -h -j $SLURM_JOB_ID -o "%o")
echo

echo "Working directory: $(pwd)"

INPUT_TEMPLATE="${input_template}"
OUTPUT_TEMPLATE="${output_json}"
INPUT_REAL="${input_real}"

# Remplacer dynamiquement le job_id dans le fichier input
sed "s/TO_BE_REPLACED_BY_SLURM/\$SLURM_JOB_ID/" "\$INPUT_TEMPLATE" > "\$INPUT_REAL"

if ! jq empty $INPUT_REAL; then
  echo "ERROR: $INPUT_REAL is not valid"
  echo FINISHED AT $(date '+%A %d %B %Y %H:%M')
  exit 1
fi

module load modules/2.3-20240529 git openblas/threaded-0.3.26
module list

echo INPUT_FILE=\$INPUT_REAL
echo OUTPUT_FILE=\$OUTPUT_JSON

echo "Permissions et contenu du dossier courant :"
ls -lh
echo

echo "Contenu du \$INPUT_REAL :"
cat "\INPUT_REAL"
echo

# Appeler ton programme principal (modifie ici selon ton code)
julia -O2 --project=/mnt/home/rbeauvallet/Many-Body-Localization  --startup-file=no --math-mode=fast /mnt/home/rbeauvallet/Many-Body-Localization/run/Cluster_running/Disorder_parallel_for_cluster.jl $INPUT_REAL
echo

echo batch file $(squeue -h -j $SLURM_JOB_ID -o "%o")
echo FINISHED AT $(date '+%A %d %B %Y %H:%M')
EOF

    chmod +x "$slurm_file"

    sbatch "$slurm_file"
done

