#!/bin/bash

## assumes htstream is available on the Pathway

start=`date +%s`
echo $HOSTNAME

inpath="00-RawData"
outpath="01-HTS_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}

for sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

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
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
