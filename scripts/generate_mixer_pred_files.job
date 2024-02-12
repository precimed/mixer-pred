#!/bin/bash
#SBATCH --job-name=MiXeR_Pred
#SBATCH --account=p697
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=16

module load R/4.0.0-foss-2020a

## Example batch command: sbatch Generate_MiXeR_Pred_files.job BIP_vs_SCZ

export MIXER_OUT_PREFIX=$1

Rscript Generate_MiXeR_Pred_files.R \
	--ref /tsd/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference/hrc_EUR_qc \
	--indir /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/MiXeR/out \
	--snps-file-prefix ${MIXER_OUT_PREFIX} \
	--outdir /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/MiXeR/MiXeR_Pred_out \
	--target-bim /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/TOP_BIPSCZvsCTRL_merged_EUR_addrs_20220127.bim \
	--plink-dir /cluster/projects/p697/users/nadinepa/Tools/plink_linux \
	--bfile /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/ref/1000G.EUR.QC.merged.all \
	--snp-thresholds 1 2 3 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300
