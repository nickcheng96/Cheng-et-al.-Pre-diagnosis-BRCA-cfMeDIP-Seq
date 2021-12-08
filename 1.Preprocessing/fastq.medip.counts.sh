#!/bin/bash
#loading modules
module load fastqc/0.11.8
module load bowtie2/2.3.5.1
module load samtools/1.9
module load picard/2.21.2
module load rstats/3.5
module load python/3.6
module load python/2.7
module load java/8
module load bedtools/2.27

#loading additional tools
export PATH=$PATH:/.mounts/labs/awadallalab/private/ncheng/softwares/BEDOPS/v2.4.38/bin
export MODULEPATH=$MODULEPATH
trimmomatic_dir="/.mounts/labs/awadallalab/private/ncheng/softwares/trimmomatic/v0.33/trimmomatic/build/bin"
bowtie2_index="/.mounts/labs/awadallalab/private/ncheng/references/hg19/Bowtie2"
bt2ref="/.mounts/labs/awadallalab/private/ncheng/references/hg19/Bowtie2/hg19_F19K16_F24B22.fa"
BEDOPS="/.mounts/labs/awadallalab/private/ncheng/softwares/BEDOPS/v2.4.38/bin"
crunch="/.mounts/labs/awadallalab/private/ncheng/softwares/ConsensusCruncher"
bwaref="/.mounts/labs/awadallalab/private/ncheng/softwares/bwa/hg19_F19K16_F24B22.fa"
bwa="/.mounts/labs/awadallalab/private/ncheng/softwares/bwa-0.7.17"
samtoolz="/.mounts/labs/awadallalab/private/ncheng/softwares/samtools/1.5/bin/samtools"
cytoband="/.mounts/labs/awadallalab/private/ncheng/softwares/ConsensusCruncher/ConsensusCruncher/hg19_cytoBand.txt"
fgbio="/.mounts/labs/awadallalab/private/ncheng/softwares/fulcrumgenomics-fgbio-7dbb0c4"
wps="/.mounts/labs/awadallalab/private/ncheng/softwares/wps_python2"
ocf="/.mounts/labs/awadallalab/private/ncheng/softwares/OCFprofiler/OCF"
my_python="/.mounts/labs/awadallalab/private/ncheng/softwares/medip_pipe/bin/python3"
picard_dir="/.mounts/labs/awadallalab/private/ncheng/softwares/picard-tools-1.119"
#location of fastq
fastq_dir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/batchall/novaseq/fastq"

#output location
outdir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/batchall/novaseq/novaseq_umitools_output_trimq"
mkdir -p $outdir
#specifying files
#runfile="AIX_0111_Ct_T_PE_61616_CM_200626_A00469_0108_BHNKKLDMXX_1_ATCCAGAG-CGACGTTA" #input file name
runfile=cur_run  
r1file=${runfile}_R1.fastq.gz
r2file=${runfile}_R2.fastq.gz
r1=$fastq_dir/$r1file
r2=$fastq_dir/$r2file
tmp=paired_trim_$r1file 
align_output=${tmp%_R1.fastq.gz}
script_dir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/batchall/scripts"
medips_script=$script_dir/medips.R
medestrand_script=$script_dir/medestrand.R
name=".R"
blist="/.mounts/labs/awadallalab/private/ncheng/softwares/ConsensusCruncher_v1/kits/IDT_duplex_sequencing_barcodes.list"
genome="hg19"
cutoff="0.7"

#start#
mkdir -p $outdir
cd $outdir

#step 0 - creating directories
mkdir -p $outdir/1_trim
mkdir -p $outdir/2_fastqc
mkdir -p $outdir/3_align_qc
mkdir -p $outdir/4_preprocess
mkdir -p $outdir/5_picard
mkdir -p $outdir/6_medips_qc/window_100 $outdir/6_medips_qc/window_200 $outdir/6_medips_qc/window_300 $outdir/6_medips_qc/window_500 $outdir/6_medips_qc/window_20 
mkdir -p $outdir/7_thalia_qc
mkdir -p $outdir/8_consensuscruncher
mkdir -p $outdir/8_qc_metrics
mkdir -p $outdir/9_medestrand/window_100 $outdir/9_medestrand/window_200 $outdir/9_medestrand/window_300 $outdir/9_medestrand/window_20 $outdir/9_medestrand/window_500
mkdir -p $outdir/10_macs2/peaks
mkdir -p $outdir/11_wps
mkdir -p $outdir/12_ocf

#step 1 - trimming adapters

#step 2 - fastqc report
fastqc $ptr1 -o 2_fastqc
fastqc $ptr2 -o 2_fastqc

#step 3 - umitools
source /.mounts/labs/awadallalab/private/ncheng/softwares/medip_venv/bin/activate


umi_tools extract -I $r1 --read2-in=$r2 --bc-pattern=NNNNN --log=processed.log --stdout=${r1}_extract --read2-out=${r2}_extract

java -jar $trimmomatic_dir/trimmomatic.jar PE \
        ${r1}_extract ${r2}_extract \
        $outdir/1_trim/paired_trim_$r1file $outdir/1_trim/unpaired_trim_$r1file $outdir/1_trim/paired_trim_$r2file $outdir/1_trim/unpaired_trim_$r2file HEADCROP:5 TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:15
       
ptr1="$outdir/1_trim/paired_trim_$r1file"
ptr2="$outdir/1_trim/paired_trim_$r2file"

#rm ${r1}_extract
#rm ${r2}_extract


bowtie2 -p 8 -x $bowtie2_index/hg19_F19K16_F24B22 \
        -1 $ptr1 \
        -2 $ptr2 \
        -S 3_align_qc/${align_output}.sam
        
samtools view -b -h 3_align_qc/${align_output}.sam  | samtools sort - -o 4_preprocess/${align_output}.sorted.bam
rm 3_align_qc/$align_output.sam

samtools index 4_preprocess/${align_output}.sorted.bam

umi_tools dedup -I 4_preprocess/${align_output}.sorted.bam --output-stats=deduplicated --output-stats=4_preprocess/${align_output}_deduplicated -S 4_preprocess/${align_output}_deduplicated.bam 


#step 4 - preprocess 
#samtools view -bS 3_align_qc/$align_output.sam | samtools sort - 4_preprocess/${align_output}.sorted
output_bam=$outdir/4_preprocess/${align_output}_deduplicated.bam
java -jar $picard_dir/MarkDuplicates.jar \
    I=$output_bam \
    O=4_preprocess/${align_output}.sorted.dedup.bam \
    M=4_preprocess/${align_output}.sorted.dedup.metrics \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=false

# step 5: get some alignment metrics
echo -e "\n~ CollectMultipleMetrics ~"
java -jar $picard_dir/CollectMultipleMetrics.jar \
    I=$output_bam \
    O=5_picard/${align_output} \
    R=$bt2ref 

echo -e "\n~ CollectGcBiasMetrics ~"
java -jar $picard_dir/CollectGcBiasMetrics.jar \
    R=$bt2ref \
    I=$output_bam \
    O=5_picard/${align_output}.gc_bias_metrics.txt \
    S=5_picard/${align_output}.summary_gc_bias_metrics.txt \
    CHART=5_picard/${align_output}.gc_bias_metrics.pdf
    
# get Thalia stats
samtools view $output_bam | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > $outdir/7_thalia_qc/${align_output}_thalia.counts
total=$(samtools view $output_bam | wc -l)
unmap=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep "^\*" | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
methyl=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep F19K16 | cut -f2); if [[ -z $methyl ]]; then methyl="0"; fi
unmeth=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep F24B22 | cut -f2); if [[ -z $unmeth ]]; then unmeth="0"; fi
pct_thalia=$(echo "scale=3; ($methyl + $unmeth)/$total * 100" | bc -l); if [[ -z $pct_thalia ]]; then pct_thalia="0"; fi
bet_thalia=$(echo "scale=3; $methyl/($methyl + $unmeth)" | bc -l); if [[ -z $bet_thalia ]]; then bet_thalia="0"; fi
echo -e "total\tunmap\tmethyl\tunmeth\tPCT_THALIANA\tTHALIANA_BETA" > $outdir/7_thalia_qc/${align_output}_thalia_summary.txt
echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_thalia\t$bet_thalia" >> $outdir/7_thalia_qc/${align_output}_thalia_summary.txt


#7 - medip windows
for window in 500 300; do
 echo -e "\n~ MEDIPS for window: $window ~"
Rscript $medips_script \
    --basedir $outdir \
    --bamfile $output_bam \
    --samplename $runfile \
    --ws $window \
    --outdir $outdir/6_medips_qc/window_${window}

 # get coverage windows
   count0=$(awk '$1 == 0' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
   count1=$(awk '$1 >= 1' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
  count10=$(awk '$1 >= 10' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
  count50=$(awk '$1 >= 50' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
 count100=$(awk '$1 >= 100' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
 echo -e "sample\tcount0\tcount1\tcount10\tcount50\tcount100" > 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt
 echo -e "$NAME\t$count0\t$count1\t$count10\t$count50\t$count100" >> 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt


 # convert to bed
 $BEDOPS/convert2bed -i wig < 6_medips_qc/window_${window}/${runfile}_medips_rpkm.wig > 6_medips_qc/window_${window}/${runfile}_${window}_medips_rpkm.bed
 $BEDOPS/convert2bed -i wig < 6_medips_qc/window_${window}/${runfile}_medips_count.wig > 6_medips_qc/window_${window}/${runfile}_${window}_medips_count.bed


 echo 'window_count' > 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt
 awk -F"\t" '$5>0' 6_medips_qc/window_${window}/${runfile}_${window}_medips_rpkm.bed | wc -l  >> 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt

 # cobble together QC data
 echo -e "\n~ Cobbling ~"
 echo -e "samples\n$align_output" > $outdir/8_qc_metrics/${runfile}.name.txt
 paste <(cat $outdir/8_qc_metrics/${runfile}.name.txt) \
       <(cut -f2-14 6_medips_qc/window_${window}/${runfile}_enrichment_data.txt) \
       <(cut -f2-3 6_medips_qc/window_${window}/${runfile}_coverage_counts.txt) \
       <(cut -f2-6 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt) \
       <(cut -f2-5 6_medips_qc/window_${window}/${runfile}_saturation_metrics.txt) \
       <(awk 'NR==7 || NR==8' 4_preprocess/${align_output}.sorted.dedup.metrics) \
       <(awk 'NR==7 || NR==8' 5_picard/${align_output}.gc_bias_metrics.txt) \
       <(awk 'NR==7 || NR==8' 5_picard/${align_output}.alignment_summary_metrics) \
       <(cut -f1 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt) \
           <(cat $outdir/7_thalia_qc/${align_output}_thalia_summary.txt) >  $outdir/8_qc_metrics/medips_${runfile}_qc_metrics_${window}.txt

done
