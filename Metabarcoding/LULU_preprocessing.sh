## This is basically what Tobias Frøslev has written in his Githut: https://github.com/tobiasgf/lulu, thanks Tobias!
# LULU preprocessing
module load blast+
makeblastdb -in ASVs.fa -parse_seqids -dbtype nucl
blastn -db ASVs.fa -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query ASVs.fa
