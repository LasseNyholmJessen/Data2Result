## Below you find a shitty written script for preprocessing data, generated by BEST/TagSteady Pipelines
## Please bear in mind that i'm not a bioinformatics person <3 
## Peace out!


### First lets demultiplex data, using AdapterRemoval. 
# Before this step please create a barcode file for each pool, using https://docs.google.com/spreadsheets/d/1R43xPdbTxKLXe5nkBNOauyHJce6YDWOuj7NT_AcwHd0/edit.
# Dont include a header in the barcode file, just copy/paste from sheet two. :)
# To use line below include paths for Fastq files and output (work)
module load AdapterRemoval
FASTQ='path/2/fastq'
WORK='path/were/you/want/your/output'
## Some times you have multiple libs/pools, therefore i run it as a loop, where each variable is a pool name with it corresponding barcode file. 
## I also normally rename raw data according to my pool name to include it in the loop
for a in 'example1' 'example2' 'example3'
  do
    AdapterRemoval --file1 $FASTQ/"$a"_1.fq.gz \
                --file2 $FASTQ/"$a"_2.fq.gz \
                --basename $WORK/"$a"/ \
                --minquality 20 --trimns --maxns 3 --trimqualities --threads 28 \
                --barcode-list $WORK/"$a".txt --barcode-mm-r1 2 --barcode-mm-r2 2
done
## clean up
for a in ''
  do
    rm $WORK/"$a"/*.settings
    rm $WORK/"$a"/*.singleton.truncated
    rm $WORK/"$a"/*.discarded
done
#Rename files
for a in 'example1' 'example2' 'example3'
  do
    echo "Yes Master, i'll rename the wierd AdapterRemoval name"
    cd "$a"
      samples=$(ls *.pair1.truncated | sed 's/.pair1.truncated//')
        for sample in ${samples}; do
          mv ${sample}.pair1.truncated ${sample}.1.fq
          mv ${sample}.pair2.truncated ${sample}.2.fq
        done
    cd ..
done

### If you don't have replicates, then skip this step. 
# Merge replicates
find *1.fq > temp
sed 's/_[1-3].1.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm temp2
sample_list=$(cat sample_list.txt)
for sample in $sample_list
  do
    cat ${sample}_1.1.fq ${sample}_2.1.fq ${sample}_3.1.fq > ${sample}.1.fq
    cat ${sample}_1.2.fq ${sample}_2.2.fq ${sample}_3.2.fq > ${sample}.2.fq
done
mkdir Single_Replicates
for sample in $sample_list
  do
    mv ${sample}_1.1.fq Single_Replicates/
    mv ${sample}_2.1.fq Single_Replicates/
    mv ${sample}_3.1.fq Single_Replicates/
    mv ${sample}_1.2.fq Single_Replicates/
    mv ${sample}_2.2.fq Single_Replicates/
    mv ${sample}_3.2.fq Single_Replicates/
done


## This is kindly borrowed from Ostaizka Aizpurua, THANKS!
# Ensure prober read orientation
module load python/v3.5.2 cutadapt/v2.6
samples=$(ls *.1.fq | sed 's/.1.fq//')
for sample in ${samples}
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq -p ${sample}_2a_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1
echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCANCAG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq -p ${sample}_2b_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq ${sample}_2b_trimmed.fq > ${sample}_1_trimmed.fq
cat ${sample}_2a_trimmed.fq ${sample}_1b_trimmed.fq > ${sample}_2_trimmed.fq
rm ${sample}_1a_trimmed.fq ${sample}_2b_trimmed.fq ${sample}_2a_trimmed.fq ${sample}_1b_trimmed.fq
done

## Continue to DADA2 pipeline
# Now data is ready to use for the DADA2 pipeline
