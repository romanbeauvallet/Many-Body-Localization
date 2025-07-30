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
    scp flatiron:/mnt/home/rbeauvallet/ceph/dataHeisenberg/ParallelDisorder/${dir}/output_h_${float_value}_beta10.json .
done
