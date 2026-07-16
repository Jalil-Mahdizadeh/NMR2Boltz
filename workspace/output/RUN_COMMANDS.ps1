$ErrorActionPreference = "Stop"
$dockerHost = "tcp://localhost:2375"
$image = "nmr2boltz:0.1.0-validated"
$hostWorkspace = (Resolve-Path .\workspace).Path
$hostRepo = (Resolve-Path .).Path

New-Item -ItemType Directory -Force .\workspace\input, .\workspace\output | Out-Null

curl.exe -L https://files.rcsb.org/download/6M6O.cif -o .\workspace\input\6M6O.cif
curl.exe -L https://files.rcsb.org/download/6M6O.pdb -o .\workspace\input\6M6O.pdb
curl.exe -L https://files.rcsb.org/download/6M6O.mr -o .\workspace\input\6M6O.mr
curl.exe -L https://files.rcsb.org/download/6M6O_nmr-data.str -o .\workspace\input\6M6O_nmr-data.str
curl.exe -L https://www.ebi.ac.uk/pdbe/entry-files/download/6m6o_validation.pdf -o .\workspace\input\6m6o_validation.pdf
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/bmr28061_3.str -o .\workspace\input\bmr28061_3.str
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/bmr28061_3.md5 -o .\workspace\input\bmr28061_3.md5
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/bmr28061.prot.fasta -o .\workspace\input\bmr28061.prot.fasta
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/ -o .\workspace\input\bmr28061_directory_index.html
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/validation/ -o .\workspace\input\bmr28061_validation_index.html
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/validation/AVS_anomalous.str.tmp -o .\workspace\input\bmr28061_AVS_anomalous.str.tmp
curl.exe -L https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr28061/validation/AVS_full.txt -o .\workspace\input\bmr28061_AVS_full.txt

docker --host $dockerHost build --tag $image .

docker --host $dockerHost run --rm --network none --read-only --tmpfs /tmp:size=64m --memory 512m --cpus 1 --user 65532:65532 --mount "type=bind,source=$hostWorkspace,target=/workspace" $image convert /workspace/input/6M6O_nmr-data.str -o /workspace/output/nmr2boltz_all --hypotheses 32 --seed 6060

docker --host $dockerHost run --rm --network none --read-only --tmpfs /tmp:size=64m --memory 512m --cpus 1 --user 65532:65532 --mount "type=bind,source=$hostWorkspace,target=/workspace" $image convert /workspace/input/6M6O_nmr-data.str -o /workspace/output/nmr2boltz_noe --origin noe --hypotheses 32 --seed 6060

docker --host $dockerHost run --rm --network none --read-only --tmpfs /tmp:size=64m --memory 1g --cpus 2 --user 65532:65532 --mount "type=bind,source=$hostWorkspace,target=/workspace" --mount "type=bind,source=$hostRepo,target=/repo,readonly" --entrypoint python $image /repo/validation/compare_ensemble.py --report /workspace/output/nmr2boltz_noe/conversion_report.json --pdb /workspace/input/6M6O.pdb --cif /workspace/input/6M6O.cif --output-dir /workspace/output/coordinate_comparison --tolerance 0.000001

docker --host $dockerHost run --rm --network none --read-only --tmpfs /tmp:size=256m --memory 1g --cpus 2 --mount "type=bind,source=$hostWorkspace,target=/workspace" --mount "type=bind,source=$hostRepo,target=/repo,readonly" --entrypoint python $image /repo/validation/stress_validation.py --iterations 100000 --seed 20260716
