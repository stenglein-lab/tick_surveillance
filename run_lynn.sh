nextflow run main.nf -c /scicomp/groups/Projects/EET_Tick_Surveillance/new_eet_scicomp.config -profile local,singularity,test \
    --outdir results_08_18_25_v206dev_targets \
    --filter_unassigned_seq Babesia,Borrelia,Borreliella,Anaplasma,Ehrlichia \
    --blast_unassigned_sequences false \
    --make_targets true \
    --primer_mix iPM-01 \
    --study_type surveillance \
    --all_refseq /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2025/02_12_25_MPAS_create_targets_process/08-13-25_MPASv206_addMakeTargets_dev/v206_dev/refseq/refseq_all.xlsx \
    --surveillance_columns /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2025/02_12_25_MPAS_create_targets_process/08-13-25_MPASv206_addMakeTargets_dev/v206_dev/refseq/surveillance_columns_surveillance_make_targets.txt \
    -work-dir /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2025/02_12_25_MPAS_create_targets_process/08-13-25_MPASv206_addMakeTargets_dev/v206_dev/work

