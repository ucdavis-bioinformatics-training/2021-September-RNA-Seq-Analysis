#!/bin/bash

## assumes htstream is available on the Path

start=`date +%s`
echo $HOSTNAME

inpath="00-RawData"
outpath="01-HTS_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}

for sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

  call="hts_Stats -L ${outpath}/${sample}/${sample}.json -N 'initial stats' \
            -1 ${inpath}/${sample}/*R1*.fastq.gz \
            -2 ${inpath}/${sample}/*R2*.fastq.gz | \
        hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'screen phix' | \
        hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'count the number of rRNA reads'\
            -r -s References/human_rrna.fasta | \
        hts_SuperDeduper -A ${outpath}/${sample}/${sample}.json -N 'remove PCR duplicates' | \
        hts_AdapterTrimmer -A ${outpath}/${sample}/${sample}.json -N 'trim adapters' | \
        hts_PolyATTrim --no-left --skip_polyT  -A ${outpath}/${sample}/${sample}.json -N 'remove polyAT tails' | \
        hts_NTrimmer -A ${outpath}/${sample}/${sample}.json -N 'remove any remaining N characters' | \
        hts_QWindowTrim -A ${outpath}/${sample}/${sample}.json -N 'quality trim the ends of reads' | \
        hts_LengthFilter -A ${outpath}/${sample}/${sample}.json -N 'remove reads < 50bp' \
            -n -m 50 | \
        hts_Stats -A ${outpath}/${sample}/${sample}.json -N 'final stats' \
            -f ${outpath}/${sample}/${sample}"

  echo $call
  eval $call
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
