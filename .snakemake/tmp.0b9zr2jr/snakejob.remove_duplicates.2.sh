#!/bin/sh
# properties = {"type": "single", "rule": "remove_duplicates", "local": false, "input": ["/project2/xinhe/alan/Cancer/testing_wgs/aligned/DO30_normal.bam"], "output": ["/project2/xinhe/alan/Cancer/testing_wgs/output/DO30_normal_sorted_dedupped.bam"], "wildcards": {"bam_name": "DO30_normal"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {}}
cd /project2/xinhe/alan/Cancer/wgs_pipeline && \
/home/selewa/.conda/envs/biotools/bin/python3.6 \
-m snakemake /project2/xinhe/alan/Cancer/testing_wgs/output/DO30_normal_sorted_dedupped.bam --snakefile /project2/xinhe/alan/Cancer/wgs_pipeline/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/xinhe/alan/Cancer/wgs_pipeline/.snakemake/tmp.0b9zr2jr /project2/xinhe/alan/Cancer/testing_wgs/aligned/DO30_normal.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /project2/xinhe/alan/Cancer/wgs_pipeline/configfile.yaml -p --nocolor \
--notemp --no-hooks --nolock --mode 2  --allowed-rules remove_duplicates  && touch "/project2/xinhe/alan/Cancer/wgs_pipeline/.snakemake/tmp.0b9zr2jr/2.jobfinished" || (touch "/project2/xinhe/alan/Cancer/wgs_pipeline/.snakemake/tmp.0b9zr2jr/2.jobfailed"; exit 1)

