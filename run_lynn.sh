nextflow run main.nf -c /scicomp/groups/Projects/EET_Tick_Surveillance/new_eet_scicomp.config -resume -profile local,singularity,test \
    --outdir results_08_18_25_v206dev_targets \
    --filter_unassigned_seq Babesia,Borrelia,Borreliella,Anaplasma,Ehrlichia \
    --blast_unassigned_sequences false \
    --make_targets true \
    --primer_mix iPM-05 \
    --study_type RESEARCH \
    --all_refseq /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2026/01-06-26_pipeline_checks_make_targets/v206_dev/refseq/refseq_all.xlsx \
    --surveillance_columns /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2026/01-06-26_pipeline_checks_make_targets/v206_dev/refseq/surveillance_columns_research_make_targets.txt \
    -work-dir /scicomp/groups/Projects/EET_Tick_Surveillance/Research/Lynn/2026/01-06-26_pipeline_checks_make_targets/v206_dev/work

