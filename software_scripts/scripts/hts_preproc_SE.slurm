#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=60
#SBATCH --mem=3000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=mrnaseq_workshop
#SBATCH --account=mrnaseq_workshop
#SBATCH --array=1-22
#SBATCH --output=slurmout/htstream_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/htstream_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail@email.com

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

inpath="00-RawData"
outpath="01-HTS_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load htstream/1.1.0

call="hts_Stats -L ${outpath}/${sample}/${sample}_htsStats.log -U ${inpath}/${sample}/*R1* | \
      hts_SeqScreener -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_SeqScreener -s References/human_rrna.fasta -r -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_AdapterTrimmer -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_PolyATTrim --no-left --skip_polyT -A ${outpath}/${sample}/${sample}.json | \
      hts_NTrimmer -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_QWindowTrim -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_LengthFilter -m 45 -A ${outpath}/${sample}/${sample}_htsStats.log | \
      hts_Stats -A ${outpath}/${sample}/${sample}_htsStats.log -f ${outpath}/${sample}/${sample}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
