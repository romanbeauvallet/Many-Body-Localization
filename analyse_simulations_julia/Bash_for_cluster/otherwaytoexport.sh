#!/bin/bash
# Variables
REMOTE_USER="rbeauvallet"
REMOTE_HOST="rustyamd2"
REMOTE_BASE_PATH="/mnt/home/rbeauvallet/ceph/dataHeisenberg/ParallelDisorder/2025_07_28"
LOCAL_DEST="../DATA_Cluster/Disorder_parallel/2025_07_28"  # Répertoire local de destination

# Créer dossier local s'il n'existe pas
mkdir -p "$LOCAL_DEST"

# Récupérer la liste des fichiers output_ à copier
ssh ${REMOTE_USER}@${REMOTE_HOST} "find ${REMOTE_BASE_PATH}/h_* -type f -name 'output_*'" > remote_output_files.txt

# Transférer chaque fichier
while read -r remote_file; do
    echo "Transferring $remote_file"
    scp "${REMOTE_USER}@${REMOTE_HOST}:${remote_file}" "$LOCAL_DEST/"
done < remote_output_files.txt
