# MiXeR-Pred
MiXeR-Pred is a novel tool for polygenic prediction that builds on an established [MiXeR](https://www.nature.com/articles/s41467-019-10310-0) framework to source genetic overlap with a secondary phenotype to inform polygenic prediction of a primary phenotype. MiXeR-Pred extends previous bivariate MiXeR models to include the functional annotations modelled in our new tool [GSA-MiXeR](https://www.medrxiv.org/content/10.1101/2022.12.08.22283159v1).

Our preprint presenting MiXeR-Pred can be found here: https://www.medrxiv.org/content/10.1101/2024.02.19.24303039v1. We ask that you please cite this paper when using our tool.

# Getting Started
## Step 1: Install GSA-MiXeR
MiXeR-Pred relies on estimates derived from our GSA-MiXeR tool. Please follow the download and installation instructions [here](https://github.com/precimed/gsa-mixer).

## Step 2: Prepare GWAS Summary Statistics
The summary statistics being used must be formated as specified [here](https://github.com/precimed/gsa-mixer?tab=readme-ov-file#input-data-formats). Most importantly the following columns must be specified:
* ``SNP`` or ``rsid`` - a marker name for each variant
* ``CHR`` - a chromosome label 
* ``BP`` or ``POS`` - genomic corrdinates for a given variant (this must be compatible with the reference build used for GSA-MiXeR (default = hg19/GRCh37))
* ``A1`` or ``EffectAllele`` - the allele used for effect estimates
* ``A2`` or ``OtherAllele``  - the non-effect/alternative allele
* ``N`` - sample size; NOTE: for case-control studies this should be the effective sample size computed as N=4/(1/ncases+1/ncontrols)
* ``Z`` - the Z-score effect estimate

NOTE: for ease of use, all input summary statistics files should be in gzip format and end in ``.sumstats.gz``. Currently, scripts are written to recognize this suffix but can be modified if needed.

## Step 3: Univariate MiXeR-Pred Run 
First, all sumstats must be processed individually to generate phenotype specific estimates. This step assumes all relevant GWAS summary statistics are in a single folder. An example script to run the univariate analysis is bellow (please modify where necessary). This script is provided in this repository at scripts/mixer_pred_univar.job and here is an example of how to batch this script for a sumstats file named ``SCZ.sumstats.gz``: ``sbatch mixer_pred_univar.job SCZ``


```
#!/bin/bash
#SBATCH --job-name=plsa2d
#SBATCH --account=p697
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=16    # remember to update --threads argument below!
#SBATCH --array=1-20
#SBATCH --chdir slurm

##################################
## The Following commands may require modification for your set-up
source /cluster/bin/jobsetup

module load singularity/3.7.1

export THREADS=16 	# note this is equal to the number of cpus used when batching above

export APPTAINER_BIND=
export SINGULARITY_BIND=
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif
md5sum ${MIXER_SIF}

export MIXER_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python         /tools/mixer/precimed/mixer.py"
export MIXER_FIGURES_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /tools/mixer/precimed/mixer_figures.py"

export ANALYSIS_ROOT=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred
export OUT_FOLDER=${ANALYSIS_ROOT}/MiXeR/out
export SUMSTATS_FOLDER=${ANALYSIS_ROOT}/sumstats

export SUMSTATS_FILE=$1

export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference
export REP=${SLURM_ARRAY_TASK_ID}
export LOADLIB_FILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/tsd_libfile_hrc_chr@.bin

export COMMON_FLAGS="${COMMON_FLAGS} --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim"
export COMMON_FLAGS="${COMMON_FLAGS} --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld"
export COMMON_FLAGS="${COMMON_FLAGS} --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE}"

export EXTRA_FLAGS_U="--exclude-ranges MHC --z1max 9.336"

## No need to modify beyond this point
##################################


# PLSA analysis (new MiXeR model)
  # remember that for PLSA analysis sumstats should be split per chromosome, as follows
${MIXER_PY} split_sumstats --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz --out ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz

  echo -e '\n=== PLSA FIT0 ===========================================================================\n'
  ${MIXER_PY} plsa --gsa-base ${COMMON_FLAGS} \
          --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
          --out ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.fit0 \
          --extract ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep${REP}.snps --hardprune-r2 0.8 \
          --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
          --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --annot-file-test ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --seed $((REP+1000)) --threads ${THREADS} ${EXTRA_FLAGS_U}

  echo -e '\n==== PLSA FIT1 ==========================================================================\n'
  ${MIXER_PY} fit1 ${COMMON_FLAGS} \
          --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz \
          --out ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.fit1 \
          --load-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.fit0.json \
          --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
          --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --extract ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep${REP}.snps --hardprune-r2 0.8 \
          --seed $((REP+1000)) --threads ${THREADS} ${EXTRA_FLAGS_U} 

  echo -e '\n==== PLSA TEST1 ==========================================================================\n'
  ${MIXER_PY} test1 ${COMMON_FLAGS} \
          --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz \
          --out ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.test1 \
          --load-params-file ${OUT_FOLDER}/${SUMSTATS_FILE}_rep${REP}.fit1.json \
          --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
          --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
          --make-snps-file \
          --seed $((REP+1000)) --threads ${THREADS} ${EXTRA_FLAGS_U} 


```

## Step 4: Bivariate MiXeR-Pred Run 
After the univariate runs are completed for your phenotypes of interest, bivariate analyses can be run for the phenotype pairs of interest. An example script to run the bivariate analysis is bellow (again, please modify where necessary). This script is provided in this repository at scripts/mixer_pred_bivar.job and here is an example of how to batch this script for the pair of sumstats files named ``BIP.sumstats.gz`` and ``SCZ.sumstats.gz``: ``sbatch mixer_pred_bivar.job BIP SCZ``

```
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

## define path to the lib file
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

```

## Step 5: Generate MiXeR-Pred Parameters for Polygenic Prediction 
After the bivariate runs are complete parameters for MiXeR-Pred must be generated. To do so you can use the ``generate_mixer_pred_files.job`` script to batch the ``generate_mixer_pred_files.R`` script. Both scripts can be found in this repository in the scripts directory. R must be pre-installed on your system (see https://cran.r-project.org/ for details). While the scripts attempt to install all required R libraries, for ease of use, you may wish to pre-insall the following libraries:  ``data.table``, ``tidyverse`` and ``argparser``. In addition, a path to a ``PLINK (v1.9)`` binary executable is needed. You can download PLINK from [here](https://www.cog-genomics.org/plink/). Note, MiXeR-Pred assumes the target sample you'd like to generate polygenic scores for has genetic data stored in PLINK binary format (.bed, .bim, .fam).

The ``generate_mixer_pred_files.job`` file passes several arguments to the R script. The follwing arguments will need to be modified:
* ``--ref`` - the path to the reference directory with the file prefix used for each chromosome
* ``--indir`` - the directory which contains all outputs from the bivariate run (Step 4)
* ``--outdir`` - the directory to place all outputs from this script
* ``--target-bim`` - the bim file for the target sample
* ``plink-dir`` - the directory containing plink binary file
* ``--bfile`` - the path and file prefix for the reference to be used in the clumping procedure (a 1000 genomes european reference file can be downloaded from [here](https://cncr.nl/research/magma/))
* ``--snp-thresholds`` - the thresholds for top number of variants (in thousands) to be used in estimating the polygenic risk scores

NOTE: the ``--snps-file-prefix`` in ``generate_mixer_pred_files.job`` does not need to be modified

Below is an example ``generate_mixer_pred_files.job`` script which can be batched using the following command given the example pair of phenotypes BIP and SCZ used previously ``sbatch generate_mixer_pred_files.job BIP_vs_SCZ``

```
#!/bin/bash
#SBATCH --job-name=MiXeR_Pred
#SBATCH --account=p697
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M 
#SBATCH --cpus-per-task=16

#load your R module
module load R/4.0.0-foss-2020a

export MIXER_OUT_PREFIX=$1

#modify paths to files/directories where necessary
Rscript Generate_MiXeR_Pred_files.R \
	--ref /tsd/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference/hrc_EUR_qc \
	--indir /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/MiXeR/out \
	--snps-file-prefix ${MIXER_OUT_PREFIX} \
	--outdir /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/MiXeR/MiXeR_Pred_out \
	--target-bim /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/TOP_BIPSCZvsCTRL_merged_EUR_addrs_20220127.bim \
	--plink-dir /cluster/projects/p697/users/nadinepa/Tools/plink_linux \
	--bfile /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/ref/1000G.EUR.QC.merged.all \
	--snp-thresholds 1 2 3 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300

```

The following outputs will be produced given the example "BIP_vs_SCZ" trait pairing:
* ``BIP_vs_SCZ_mixer_estimates.gz`` - a file with MiXeR-Pred estimates for all variants in the reference as well as both input summary statistics
* ``BIP_vs_SCZ_mixer_estimates_target_filtered.gz`` - a file with MiXeR-Pred estimates for all ``BIP_vs_SCZ_mixer_estimates.gz`` file also found in the target sampe .bim file
* ``BIP_vs_SCZ_mixer_target_filtered.clumped`` - independent significant variants used for polygenic score analysis
* ``BIP_vs_SCZ_mixer_target_filtered.log`` and BIP_vs_SCZ_mixer_target_filtered.nosex - standard output files produced from the plink clumping step.
* ``BIP_vs_SCZ_clump_topSNPs`` - a directory with all files used to generate polygenic risk scores at the pre-selected number of top variants
  * these files are named ``BIP_vs_SCZ_top_(n)K.mixer.gz``, where "(n)" represents the number of included top independent significant variants in the thousands.

All the gzipped output files will have the following columns:
* ``RSID``- a marker name for each variant
* ``CHR`` - a chromosome label 
* ``POS`` - genomic corrdinates for a given variant (in build hg19/GRCh37)
* ``A1`` - the effect allele
* ``A2`` - the non-effect/alternative allele
* ``ED1`` - expected delta for the primary phenotype given the secondaty phenotype (i.e., weight used for the the first listed phenotype (BIP))
* ``ED2`` - expected delta for the secondary phenotype given the primary phenotype (i.e., weight used for the the second listed phenotype (SCZ))
* ``expED1`` - exponentiated negative of ED1 (i.e., thresholding value used for the the first listed phenotype (BIP))
* ``expED2`` - exponentiated negative of ED2 (i.e., thresholding value used for the the second listed phenotype (SCZ))

NOTE: polygenic scores are generated using ED1 (or ED2) as weights and expED1 (or expED2) as thresholds similar to using beta values as weights and p-values as thresholds. This means weights and thresholds are generated for both phenotypes in a given pair of traits.

## Step 6: Generate Polygenic Scores
The output files found in the "*_clump_topSNPs" directory can be used to generate polygenic scores. In the MiXeR-Pred manuscript (link coming soon), we use PRSice2 to generate polygenic scores. PRSice2 can be downloaded [here](https://choishingwan.github.io/PRSice/). Alternatively plink can be used to generate polygenic scores using the ``--score`` command (https://www.cog-genomics.org/plink/1.9/score).

An example script to generate polygenic scores at multiple thresholds using PRSice2 is below. NOTE: ``--stat ED1`` ``--pvalue expED1`` sets the MiXeR-Pred weights and thresholds, respectively.
```
#!/bin/bash
#SBATCH --job-name=Pred
#SBATCH --account=p697_norment_dev
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000M

# load R module
module load R/4.0.0-foss-2020a

# Top variant thresholds to use for PRS
declare -a Nmax_arr2=(1 2 5 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300)
for Nmax in "${Nmax_arr2[@]}"
do
# define base files - output files from Generate_MiXeR_Pred_files.job script
INP_BASE=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/MiXeR/MiXeR_Pred_out/BIP_vs_SCZ_clump_topSNPs/BIP_vs_SCZ_top_${Nmax}K.mixer.gz

  # Path to phenotype file used for PRSice2
	INP_PHENO=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/TOP_psych_DIAGnSYMP_BIPvsCTRL_20220127.pheno
  # Pathe to covariate file used for PRSice2
	INP_COV=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/TOP_psych_DIAGnSYMP_20211115.cov
  # Path to target sample plink files with prefix for binary files (.bed, .bim. .fam files)
	INP_TARGET=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/TOP_BIPSCZvsCTRL_merged_EUR
	#Output file prefix
	flag=/cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/mixer_prs/BIP_vs_SCZ_${Nmax}K
						
		#===ANALYZING===
		if [ -f "$INP_BASE" ]; 
		then
      # Standard commands to run PRSice2
			Rscript /cluster/projects/p697/users/nadinepa/Tools/PRSice_linux/PRSice.R \
				--prsice /cluster/projects/p697/users/nadinepa/Tools/PRSice_linux/PRSice_linux \
				--dir /cluster/projects/p697/users/nadinepa/MiXeR_Pred/FINAL_MiXeR_Pred/PRS/mixer_prs \
				--all-score \
				--model add \
				--no-clump \
				--bar-levels 1 --fastscore \
				--base $INP_BASE\
				--snp RSID --A1 A1 --A2 A2 \
				--stat ED1 --pvalue expED1 --beta \
				--target $INP_TARGET \
				--pheno-file $INP_PHENO \
				--pheno-col "psych"\
				--binary-target "T" \
				--cov-file $INP_COV \
				--cov-col sex,genetics_batch,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
				--cov-factor sex,genetics_batch \
				--out ${flag}
		else
			echo "Your Base $INP_BASE does not exist!\n" > ${flag}.err
		fi
done

```

