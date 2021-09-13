#!/bin/bash

#SBATCH --job-name=salmon_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=60
#SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/salmon-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/salmon-index_%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail@email.com

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}
cd ${outpath}

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.pc_transcripts.fa.gz
gunzip gencode.vM25.pc_transcripts.fa.gz
PC_FASTA="gencode.vM25.pc_transcripts.fa"
INDEX="salmon_gencode.vM25.index"

module load salmon
call="salmon index -i ${INDEX} -k 31 --gencode -p 8 -t ${PC_FASTA}"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
