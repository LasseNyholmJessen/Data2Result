module load AdapterRemoval
FASTA='/groups/hologenomics/jagerbo/data/5-HappyFish/2-16S_profiling/0-data'
WORK='/groups/hologenomics/jagerbo/data/5-HappyFish/2-16S_profiling/1-Demux'
for a in 'D9' 'D14'
  do
    AdapterRemoval --file1 $FASTA/"$a"_1.fq.gz \
                --file2 $FASTA/"$a"_2.fq.gz \
                --basename $WORK/"$a"/ \
                --minquality 20 --trimns --maxns 3 --trimqualities --threads 28 \
                --barcode-list $WORK/"$a".txt --barcode-mm-r1 2 --barcode-mm-r2 2
done

## clean up
for a in 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D17' 'D18'
  do
    rm $WORK/"$a"/*.settings
    rm $WORK/"$a"/*.singleton.truncated
    rm $WORK/"$a"/*.discarded
done
#Rename files
for a in 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D17' 'D18'
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

for a in 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D17' 'D18'
  do
    echo "Yes Master, i'll move samples from $a to the new directory"
    cp $WORK/"$a"/*.fq $WORK/All/
done

# Merge replicates
cd All
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

# Clean up
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
