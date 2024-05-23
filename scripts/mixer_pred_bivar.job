#!/bin/bash
#SBATCH --job-name=bivar
#SBATCH --account=p697
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=16    # remember to update --threads argument below!
#SBATCH --array=1-20
#SBATCH --chdir slurm

##################################
## The Following commands may require modification for your set-up
## NOTE: These commands can be identical to the univariate analysis

source /cluster/bin/jobsetup

# modify the number of threads to match the number of cpus used for batching above
export THREADS=16 

# Define a working directory as ANYSIS_ROOT and ouput directory as OUT_FOLDER
export ANALYSIS_ROOT=/cluster/projects/p697/users/nadinepa/MiXeR_Pred
export OUT_FOLDER=${ANALYSIS_ROOT}/MiXeR/out

# Provide the location of the summary statistics directory
export SUMSTATS_FOLDER=${ANALYSIS_ROOT}/sumstats

## load singularity module for batching
module load singularity/3.7.1 

## provide the path to the GSA-MiXeR singularity container
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif

# define MIXER_PY command which executes mixer.py script within the singularity container
export MIXER_FIGURES_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /tools/mixer/precimed/mixer_figures.py"

# deifne the path to the ld reference directory
export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference

export APPTAINER_BIND=
export SINGULARITY_BIND=
export REP=${SLURM_ARRAY_TASK_ID}

## NEED ALEX'S HELP!! latest GSA-MIXER no longer defines this argument
export LOADLIB_FILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/tsd_libfile_hrc_chr@.bin

## No need to modify beyond this point
##################################

export SUMSTATS_FILE=$1
export SUMSTATS2_FILE=$2

export COMMON_FLAGS="${COMMON_FLAGS} --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim"
export COMMON_FLAGS="${COMMON_FLAGS} --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld"
export COMMON_FLAGS="${COMMON_FLAGS} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE}"

export EXTRA_FLAGS_B="--exclude-ranges MHC --z1max 9.336 --z2max 9.336"

# PLSA analysis (new MiXeR model)
  echo -e '\n==== PLSA FIT2 ==========================================================================\n'
  ${MIXER_PY} fit2 ${COMMON_FLAGS} \
          --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz \
          --trait2-file ${SUMSTATS_FOLDER}/${SUMSTATS2_FILE}.sumstats.gz \
          --trait1-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.fit1.json \
          --trait2-params-file ${OUT_FOLDER}/${SUMSTATS2_FILE}_rep${REP}.fit1.json \
          --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
          --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --make-snps-file \
          --out ${OUT_FOLDER}/${SUMSTATS_FILE}_vs_${SUMSTATS2_FILE}_rep${REP}.fit2 \
          --extract ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep${REP}.snps --hardprune-r2 0.8 \
          --seed $((REP+1000)) --threads ${THREADS} ${EXTRA_FLAGS_B}

  echo -e '\n==== PLSA TEST2 ==========================================================================\n'
  ${MIXER_PY} test2 ${COMMON_FLAGS} \
          --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz \
          --trait2-file ${SUMSTATS_FOLDER}/${SUMSTATS2_FILE}.sumstats.gz \
          --load-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_vs_${SUMSTATS2_FILE}_rep${REP}.fit2.json \
          --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
          --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --make-snps-file \
          --out ${OUT_FOLDER}/${SUMSTATS_FILE}_vs_${SUMSTATS2_FILE}_rep${REP}.test2 \
          --seed $((REP+1000)) --threads ${THREADS} ${EXTRA_FLAGS_B} 
