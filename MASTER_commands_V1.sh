


#########################
# Read trimming 
#########################

#trassembly_phylogenomics.sh
fastp --in1 GK1_R1.fastq --in2 GK1_R2.fastq --out1 GK1_R1_trimmed.fastq.gz --out2 GK1_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT



#trassembly_auto_attempt.sh
fastp --in1 D441_3_08144AAE_GCCGCAAGTG-TTACCGCAAT_R1_001.fastq.gz --in2 D441_3_08144AAE_GCCGCAAGTG-TTACCGCAAT_R2_001.fastq.gz --out1 D44_R1_trimmed.fastq.gz --out2 D44_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 10D*_R1_001.fastq.gz --in2 10D*_R2_001.fastq.gz --out1 10D_R1_trimmed.fastq.gz --out2 10D_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 13G*_R1_001.fastq.gz --in2 13G*_R2_001.fastq.gz --out1 13G_R1_trimmed.fastq.gz --out2 13G_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 16Ckin*_R1_001.fastq.gz --in2 16Ckin*_R2_001.fastq.gz --out1 16Ckin_R1_trimmed.fastq.gz --out2 16Ckin_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT



#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do fastp --in1 "$i"*_R1_001.fastq.gz --in2 "$i"*_R2_001.fastq.gz --out1 "$i"_R1_trimmed.fastq.gz --out2 "$i"_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ; done


#########################
# ASSEMBLY 
#########################
# (originating from different batches, reads from different cultures were assembles at different times, although with technically similar commands. For the record, I keep all of them as they were ran, and the preceeding line corresponds to the name of the script which contains original commands ran on a given batch)

#trassembly_phylogenomics.sh
rnaspades -1 GK1_R1_trimmed.fastq.gz -2 GK1_R2_trimmed.fastq.gz --ss-fr -o /scratch/data2/dzavadska/RNAspades/Gemkin/
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file 


#alt_assemblers_trassembly_auto_attempt.sh - original source of the following line
for i in $(echo "D44 10D 13G 16Ckin") ; do nohup rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done


#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done





#########################
# DECONTAMINATION 
#########################

##############
# PhyloFlash #
##############

#trassembly_phylogenomics.sh
/home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib Gemkin -read1 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R1_trimmed.fastq.gz -read2 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R2_trimmed.fastq.gz -almosteverything -readlimit 10000000 -readlength 86

## for "D44 10D 13G 16Ckin" , PhyloFlash was ran by Daniel

#trassembly_Jun2024_auto_attempt.sh - original source of the following li


#########################
# Read trimming 
#########################

#trassembly_phylogenomics.sh
fastp --in1 GK1_R1.fastq --in2 GK1_R2.fastq --out1 GK1_R1_trimmed.fastq.gz --out2 GK1_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT



#trassembly_auto_attempt.sh
fastp --in1 D441_3_08144AAE_GCCGCAAGTG-TTACCGCAAT_R1_001.fastq.gz --in2 D441_3_08144AAE_GCCGCAAGTG-TTACCGCAAT_R2_001.fastq.gz --out1 D44_R1_trimmed.fastq.gz --out2 D44_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 10D*_R1_001.fastq.gz --in2 10D*_R2_001.fastq.gz --out1 10D_R1_trimmed.fastq.gz --out2 10D_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 13G*_R1_001.fastq.gz --in2 13G*_R2_001.fastq.gz --out1 13G_R1_trimmed.fastq.gz --out2 13G_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

fastp --in1 16Ckin*_R1_001.fastq.gz --in2 16Ckin*_R2_001.fastq.gz --out1 16Ckin_R1_trimmed.fastq.gz --out2 16Ckin_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT



#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do fastp --in1 "$i"*_R1_001.fastq.gz --in2 "$i"*_R2_001.fastq.gz --out1 "$i"_R1_trimmed.fastq.gz --out2 "$i"_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ; done


#########################
# ASSEMBLY 
#########################
# (originating from different batches, reads from different cultures were assembles at different times, although with technically similar commands. For the record, I keep all of them as they were ran, and the preceeding line corresponds to the name of the script which contains original commands ran on a given batch)

#trassembly_phylogenomics.sh
rnaspades -1 GK1_R1_trimmed.fastq.gz -2 GK1_R2_trimmed.fastq.gz --ss-fr -o /scratch/data2/dzavadska/RNAspades/Gemkin/
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file 


#alt_assemblers_trassembly_auto_attempt.sh - original source of the following line
for i in $(echo "D44 10D 13G 16Ckin") ; do nohup rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done


#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done





#########################
# DECONTAMINATION 
#########################

##############
# PhyloFlash #
##############

#trassembly_phylogenomics.sh
/home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib Gemkin -read1 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R1_trimmed.fastq.gz -read2 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R2_trimmed.fastq.gz -almosteverything -readlimit 10000000 -readlength 86

## for "D44 10D 13G 16Ckin" , PhyloFlash was ran by Daniel

#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do /home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib "$i" -read1 "$i"*_R1_001.fastq.gz -read2 "$i"*_R2_001.fastq.gz -almosteverything -readlimit 10000000 -readlength 150 ; done




#copy phyloflash outputs

# cp /scratch/data1/agalvez/daryna/phyloflash/G65133.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data1/agalvez/daryna/phyloflash/RhMon.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data2/dzavadska/tr_ass_blastdbs/*.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash

# save contaminant genome IDs to cont_genomes.txt

#Download all putative phyloflash contaminant genomes, unzip them and get genomic fnas

cat cont_genomes.txt | while read i ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done 
#for i in $(echo "GCA_030717845.1" "GCA_003259155.1") ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done

#copy all assemblies to denctm folder and create blast databases

cp ./* ../decontam/ass_blastdbs 
cd ../decontam/ass_blastdbs 
for file in *.fasta; do makeblastdb -in $file -dbtype nucl; done


#2) Now BLAST!!!!!!

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query queries.fasta -db /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/ass_blastdbs/"$i"_transcripts.fasta -outfmt 6 -out "$i"_blast_ssu_vs_ass.txt ; done

#manually select putative true matches from the list and create lists of subject contigs to be retrieved and subsequently blasted vs NCBI

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -r -f "$i"_retrieve_cont.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_cont_to_blast.fasta ; done

# the commented versions below hit ncbi timeout; used web_blast.pl by Daniel, which implements megablast
#for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query "$i"_cont_to_blast.fasta -db nt -remote -outfmt "6 qseqid sseqid pident length evalue stitle sskingdoms sscinames staxids" -max_target_seqs 1 -out "$i"_blast_contigs_vs_ncbi.txt ; done

#/home/dzavadska/ncbi-blast-2.16.0+/bin/blastn -query D44_cont_to_blast.fasta -db /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -outfmt "6 qseqid sseqid pident length evalue stitle" -max_target_seqs 1 -out D44_blast_contigs_vs_ncbi.txt

#blastall -p blastn -d /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -i test.fasta -o testout.txt


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do perl ./web_blast.pl megablast nt Tabular 1 0 0 "$i"_cont_to_blast.fasta> "$i"_cont_to_blast.fasta ; done
# separately for G65133 because otherwise exceeds memory limit perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_1.fasta > G65133_blast_contigs_vs_ncbi_1.txt
# perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_2.fasta > G65133_blast_contigs_vs_ncbi_2.txt
for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do print "$i" ; grep "Query=" "$i"_blast_contigs_vs_ncbi.txt ; grep "Value" -A 2 "$i"_blast_contigs_vs_ncbi.txt | grep "%" ; done
# grep "Query=" G65133_blast_contigs_vs_ncbi*
# grep "Value" -A 2 G65133_blast_contigs_vs_ncbi* | grep "%"

## manually copy all top hit ids and svae into ncbi_hit_ids.txt

cat ncbi_hit_ids.txt | while read i ; do esummary -db nuccore -id "$i" | xtract -pattern DocumentSummary -element TaxId >> TAXID_ncbi_hit_ids.txt ; done
sort -u TAXID_ncbi_hit_ids.txt > TAXID_ncbi_hit_ids_uniq.txt


cat TAXID_ncbi_hit_ids_uniq.txt | xargs -n 1 -I{} sh -c 'sleep 1; esearch -db assembly -query "{}[TAXID]" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession | head -n 5 ' >> ncbi_blast_cont_genomes.txt

mv ncbi_blast_cont_genomes.txt ../contaminant_genomes/

#in ../contaminant_genomes/
cat ncbi_blast_cont_genomes.txt cont_genomes.txt | uniq >> all_cont_genomes_uniq.txt
###Download all contaminant genomes

cat all_cont_genomes_uniq.txt | while read i ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done 

#make blast databases from the downloaded genomes
for i in ./*.fna ; do makeblastdb -in "$i" -dbtype nucl -out "$i" ; done

#BLAST transcriptomes vs contaminant genome databases
/home/drichter/bin/perl/run_BLAST.pl -query /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/ass_blastdbs -database /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/contaminant_genomes -output /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/transcr_vs_cont_genome_out


#subset 60000 random hits
cat ./transcr_vs_cont_genome_out/*decontam_refG_out6.txt > all_decontam_refG.csv
shuf -n 60000 all_decontam_refG.csv > subset_all.csv

#produce a histogram of hits in R, due to huge number of hits, just take 60000 random sample of all hits and plot those using the script in R
# /home/dzavadska/Data/COMPARATIVE_DATASET/phylogen_DECONTAM/blast_res_viz.Rmd


##Remove putative contaminant contigs
for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -v -f "$i"_contaminant_contig_list.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_transcripts_clean.fasta ; done

#double-checking number of contigs before and after decontam

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do echo $i ; grep ">" "$i"_transcripts_clean.fasta | wc -l ; grep ">" ../ass_blastdbs/"$i"_transcripts.fasta | wc -l ; done

# Gemkin 31846 32513
# RhMon 60289 62845
# G65133 209809 248992
# D44 33506 34457
# 16Ckin 33503 33787
# 13G 73868 94536
# 10D 65003 66351

#wcleaner

#wcleaner settings that worked ages ago /scratch/data2/dzavadska/wcleaner_input/settings.yml 

# Batch 1
cp /scratch/data2/dzavadska/wcleaner_input/summer_2023_batch/BEAP* ./

/scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/clean_from_contamination_trass/*_clean.fasta


for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_R1_trimmed.fastq "$i"trimmed_1.fastq ; done
for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_R2_trimmed.fastq "$i"trimmed_2.fastq ; done
for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_transcripts_clean.fasta "$i"trimmed.fasta ; done


for i in $(echo "D44" "16Ckin" "13G" "10D") ; do sed -i 's/_length.*//g' "$i"trimmed.fasta ; done


#Batch 2
fastp --in1 MP2_06939AAF_GCATAGCTAC-AATCGAGGCC_R1_001.fastq.gz --in2 MP2_06939AAF_GCATAGCTAC-AATCGAGGCC_R2_001.fastq.gz --out1 MP2_R1_trimmed.fastq.gz --out2 MP2_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

sudo cp /home/mskamnelou/data/MAOP2_beap0263/MP2_R*_trimmed.fastq.gz

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/fastp/trimmed_D57RNA_06935AAF_AAGAGGATCT-GCTCCTGACG_R*_001.fastq ./
mv trimmed_D57RNA_06935AAF_AAGAGGATCT-GCTCCTGACG_R1_001.fastq D57_R1_trimmed.fastq

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/fastp/trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R*_001.fastq ./
mv trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R1_001.fastq G21_R1_trimmed.fastq
mv trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R2_001.fastq G21_R2_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/ChY_BEOP0192_Choano_yana_until_trinity/3_Fastp/trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq.gz ./
gunzip trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq.gz
mv trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq ChY_R1_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/ChY_BEOP0192_Choano_yana_until_trinity/3_Fastp/trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R2_001.fastq ./
mv trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R2_001.fastq ChY_R2_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/FRCH_BEAPO279_FRESCHO-3_until_trinity/3_Fastp/trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R*_001.fastq ./
mv trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R1_001.fastq FRCH_R1_trimmed.fastq
mv trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R2_001.fastq FRCH_R2_trimmed.fastq

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/Chy_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/d57_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/FRCH_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/g21_1st_decont_trinity_out_dir.Trinity.fasta ./


for i in $(echo "RhMon" "G65133" "MP2") ; do sed -i 's/_length.*//g' "$i"trimmed.fasta ; done
for i in $(echo "G21" "ChY" "FRCH" "D57") ; do sed -i 's/ len.*//g' "$i"trimmed.fasta ; done


###############################
# CROSS-CONTAMINATION REMOVAL #
###############################
#RUNNNING WCLEANER



source /home/drichter/lib/external/WinstonCleaner/venv/bin/activate


#Run data prep; takes some time ~o/n
winston_prepare_data
#Run cleaner; takes 0 time
winston_find_contaminations

#deactivate wcleaner 
deactivate

#WCLEANER OUT

#Check the number of deleted contigs from *deleted.fasta files

for file in *deleted.fasta ; do echo $file ; grep '>' $file | wc -l ; done

# 10Dtrimmed_deleted.fasta
# 325
# 13Gtrimmed_deleted.fasta
# 915
# 16Cinittrimmed_deleted.fasta
# 56728
# 16Ckintrimmed_deleted.fasta
# 712
# D44trimmed_deleted.fasta
# 166

# G65133trimmed_deleted.fasta
# 146
# RhMontrimmed_deleted.fasta
# 2771


#check distribution of contig length (after decontam) - just in case

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' "$i"_doubleclean_contigs.fasta > "$i"_len_distr.txt ; done





#########################
# PROTEIN PREDICTION 
#########################



###########################
#Running TransDecoder


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do mv "$i"*.fasta "$i"_doubleclean_contigs.fasta ; done
cp Gemkin_doubleclean_contigs.fasta Gemkin_doubleclean_contigs_RC.fasta

for file in ./Gemkin*_doubleclean_contigs_RC.fasta ; do sudo TransDecoder.LongOrfs -t $file ; done
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file - therefor, ThransDecoder ran in a different mode

for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.LongOrfs -t $file -S ; done


for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.Predict -t $file ; done

#########################
# BUSCO COMPLETENESS 
#########################

#running BUSCO

# cp ./*transdecoder.pep ./busco
# cp ./*transdecoder.cds ./busco

for file in *transdecoder.pep ; do busco -i $file -l eukaryota_odb10 -o "$file"pep -m proteins --cpu 8 ; done
for file in *transdecoder.pep ; do busco -i $file -l euglenozoa_odb10 -o "$file"EUGLENOZOApep -m proteins --cpu 8 ; done
















####################################################################
####################################################################
####################################################################
# Orthology inference and verification
####################################################################
####################################################################
####################################################################


# OrthoFinder 45 taxa, CD-hit 90%, diamond_ultra_sens , msa

/home/agalvez/programs/OrthoFinder/orthofinder -f NR_proteins_90_45taxa -S diamond_ultra_sens -M msa -n "Results_diamond_ultra_sens_FASTTREE" -a 16 -t 16


#to obtain statistics per species 

# cp ./SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv  ./Summaries_SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa_Results_diamond_ultra_sens_FASTTREE_Statistics_PerSpecies.tsv

# cp ./SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Orthogroups/Orthogroups.GeneCount.tsv ./Summaries_SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa_Results_diamond_ultra_sens_FASTTREE_Orthogroups.GeneCount.tsv




## to produce "OG_summary_NR_proteins_45_taxa_diamond_ultra_sens.csv" , run Orthofinder_stats_UNIVERSAL.Rmd

#filtering in R from OG_summary:

#OG_local <- fread("OG_summary_NR_proteins_45_taxa_diamond_ultra_sens.csv")
#subset <- OG_local[which(OG_local$mean<=1 & OG_local$sd<=1 & OG_local$sp_repr>25),]
#double_subset <- subset[which(subset$mean>=0.9 & subset$sd<=0.5 & subset$sp_repr>35),]
#fwrite(double_subset, "double_subset_NR_proteins_45_taxa_diamond_ultra_sens.csv")

## 494 OG names saved into OG_double_subset.txt


#cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing
cat OG_double_subset.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Gene_Trees/"$i"* ./ ; done

cat OG_double_subset.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Orthogroup_Sequences/"$i"* ./ ; done

# renaming sequence headers

cat OG_double_subset.txt | while read i ; do sed 's/;/_/g' ./"$i".fa > ./"$i"_charsubs.fa ; done
cat OG_double_subset.txt | while read i ; do phykit tip_labels ./"$i"_tree.txt > ./"$i"_tiplabels.txt ; done
cat OG_double_subset.txt | while read i ; do grep ">" ./"$i"_charsubs.fa | sed 's/>//g' > ./"$i"_fastalabels.txt ; done
cat OG_double_subset.txt | while read e ; do while read i; do sed -i "s@"$i"@$(grep -F "$i" ./"$e"_tiplabels.txt)@g" ./"$e"_charsubs.fa ; done < "$e"_fastalabels.txt ; done
cat OG_double_subset.txt | while read e ; do mv ./"$e"_charsubs.fa ./"$e"_sp_renamed.fasta ; done


#aligning, trimming and building tree for SG-tree check
cat OG_double_subset.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 16 --reorder ""$i"_sp_renamed.fasta"  > ""$i"_sp_renamed_MAFFT.fasta" ; done

#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/TRIMMED_ASSESSED
cat ./OG_double_subset.txt | while read i ; do trimal -in "$i"_sp_renamed_MAFFT.fasta -out "$i"_mafft_gappyout.fasta -gappyout -fasta ; done

#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/IQTREE_sg/gappyout_iqtree/
for i in ./*_mafft_gappyout.fasta ; do iqtree -bb 1000 -seed 1234 -nt 24 -s "$i" ; done
# generates the set of trees used FOR THE ACTUAL CHECKS; *.contree


##### After the manual check

#copy trees and alignments that were checked from dazhbog  
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/downstream_orthogroups/orthosnap/Genes_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/downstream_orthogroups/orthosnap/Trees_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN
#copy the initial OG sequences to the new directory
# cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/*_sp_renamed.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/TO_CHECK_UPDATED_no_orthoSNAP/

#replacing | in fasta headers with _
for i in ./TO_CHECK_UPDATED_no_orthoSNAP/* ; do sed -i 's/|/_/g' $i ; done

#sed 's/|/_/g' ./removing_genes/OG0003302_sp_renamed.fasta #debug, check Sspe gene names
#check 
grep -F "|" ./TO_CHECK_UPDATED_no_orthoSNAP/*


# deleting genes

cat Genes_to_delete.tsv | cut -f 1 | uniq > uniq_OG_list.txt

cat uniq_OG_list.txt | while read i ; do grep "$i" Genes_to_delete.tsv | cut -f 2 > ./removing_genes/"$i"_to_delete.txt ; done

cat uniq_OG_list.txt | while read i ; do seqkit grep -n -v -f ./removing_genes/"$i"_to_delete.txt ./TO_CHECK_UPDATED_no_orthoSNAP/"$i"_sp_renamed.fasta -o "$i"_AFTERCHECK.fasta ; done


## Check if everything was removed
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep "$e" "$i"_AFTERCHECK.fasta ; echo "$i" ; done ; done
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep -F "$e" "$i"_AFTERCHECK.fasta ; done ; done

#manually remove few sequences that failed to be removed by seqkit for whatever reason
#>RhMon_Gene.33090__NODE_18304__g.33090__m.33090
#OG0002451
#>Cruzella_marina_Gene.84135__NODE_19292_length_2400_cov_477.822643_g7898_i0__g.84135__m.84135
#OG0002445
#>G65133_Gene.62375__NODE_45689__g.62375__m.62375
#OG0002452
#>RhMon_Gene.25779__NODE_13194__g.25779__m.25779
#OG0002452


# seqkit grep -v -f ./removing_genes/OG0003302_to_delete.txt ./TO_CHECK_UPDATED_no_orthoSNAP/OG0003302_sp_renamed.fasta

########## REALIGNING

mv ./*_AFTERCHECK.fasta ./realign/
cd ./realign/
cat ../uniq_OG_list.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 100 --reorder ""$i"_AFTERCHECK.fasta"  > ""$i"_AFTERCHECK_realigned.fasta" ; done



cat ./uniq_OG_list.txt | while read i ; do cp ./realign/"$i"_AFTERCHECK_realigned.fasta ./to_modeltest/ ; done

## removing trees marked for complete removal
cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest
cat ../Trees_to_delete.tsv | while read i ; do rm ./"$i"_AFTERCHECK_realigned.fasta ; done


#trimming the re-aligned ones

#for i in ./*AFTERCHECK_realigned.fasta  ; do trimal -in "$i" -out "$i"_gappyout.fasta -gappyout -fasta ; done
#WARNING: Removing sequence 'Gene.47917__CAMNT_0016070987__g.47917__m.47917' composed only by gaps
#WARNING: Removing sequence 'Gene.34531__NODE_19382__g.34531__m.34531' composed only by gaps

cat ./Genes_to_delete.tsv | cut -f 1 | sort -u | > OGs_realigned.txt

##adding alignments which are meant to be kept as they were initially
cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/*_sp_renamed_MAFFT.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/tmp/
cat ./Trees_to_delete.tsv | while read i ; do rm ./tmp/"$i"* ; done
cat ./uniq_OG_list.txt | while read i ; do rm ./tmp/"$i"* ; done
cp ./tmp/* ./to_modeltest/



#
ls | sed 's/_.*//g' > ../All_OGs_kept.txt
#
cat ../All_OGs_kept.txt | while read i ; do trimal -in ./"$i"* -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done

#WARNING: Removing sequence 'G65133_Gene.88879__NODE_82234__g.88879__m.88879' composed only by gaps

#renaming for concatenation

#cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/FASTTREE_draft/gappyout_fasttree/species_names.txt /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/
for file in ./*_gappyout.fasta ; do cat ../species_names.txt | while read i ; do sed -i "s/${i}_/${i}|/g" "$file" ; done ; done

## checking if there are duplicates
for i in ./*_gappyout.fasta ; do echo $(grep ">" "$i" | sed 's/|.*//g' | sort | uniq -cd ) ; echo "$i" ; done

## checking if every header got substituted with | species separator
for i in ./*_gappyout.fasta ; do echo "$i $(grep -c '^>' "$i") $(grep -c '|' "$i")" ; done

#re-doing a few duplicates detected, just manually on dazhbog
#2 >G65133
#./OG0002646_toMT_gappyout.fasta
#2 >Strigomonas_culicis
#./OG0002823_toMT_gappyout.fasta
#2 >G65133
#./OG0002831_toMT_gappyout.fasta

scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/OG0002831* ./
for i in ./* ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 100 --reorder "$i"  > ""$i"_AFTERCHECK_realigned.fasta" ; done
for i in $(echo "OG0002646" "OG0002823" "OG0002831")  ; do /home/dzavadska/Data/soft/trimAl/source/trimal -in ./"$i"_AFTERCHECK_realigned.fasta_AFTERCHECK_realigned.fasta -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done
scp ./*gappyout.fasta nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/



# for i in $(echo "OG0002650" "OG0003130" "OG0003260" "OG0003298")  ; do trimal -in ./"$i"*_realigned.fasta -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done
for file in $(echo  "OG0002646" "OG0002823" "OG0002831") ; do cat ../species_names.txt | while read i ; do sed -i "s/${i}_/${i}|/g" "$file"*_gappyout.fasta ; done ; done

#substitute specifically >g31 to Willaertia_magna|g31 - for some reason I dont really remeber, this is some strange exception.
grep ">g31" ./*_gappyout.fasta
#./OG0002815_toMT_gappyout.fasta:>g31
sed -i 's/>g31/>Willaertia_magna|g31/g' ./OG0002815_toMT_gappyout.fasta
grep ">" ./OG0002815_toMT_gappyout.fasta

#####################################################################################################
############################# CONCATENATING ######################################
#####################################################################################################
ls *_gappyout.fasta > alignments.txt

#replacing | with " "
for i in *_gappyout.fasta ; do cp "$i" ./ "$i"_BACKUP.fasta ; done
# debug
#sed 's/|/ /g' OG0003327_toMT_gappyout.fasta_BACKUP.fasta | head
for i in *_gappyout.fasta ; do sed -i 's/|/ /g' "$i" ; done

pk_create_concat -a alignments.txt -p CONCAT

#check if concatenation was done correctly; there should be only the names of taxa, same number as that of intital taxa in a dataset
grep ">" CONCAT.fa

#Partition file output: CONCAT.partition
#Concatenated fasta output: CONCAT.fa
#Occupancy report: CONCAT.occupancy

cp ./CONCAT* ../../../CONCAT_AND_GIANT_TREES/model_selection/




########################################################################################################
# TREE RUNS

#######################################
######### PHYLOBAYES SBATCH ###########

# As for choosing the right parallelization scheme, the most straightforward approach is to set n equal to the number of cores of a given node.

      
#######################################
######### PHYLOBAYES SBATCH ###########

#!/bin/bash
#----------------------------------------------------
# PHYLOBAYES SBATCH based on CESGA example 
#----------------------------------------------------
#SBATCH -J cesgaPHYTESTconcat_pb_catgtr_ch2            # Job name
#SBATCH -o cesgaPHYTESTconcat_pb_catgtr_ch2_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e cesgaPHYTESTconcat_pb_catgtr_ch2_%j.e       # Name of stderr output file(%j expands to jobId)
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB           # Cores per task requested
#SBATCH -t 124:00:00         # Run time (hh:mm:ss) - 2 days     
      
module load cesga/2020 gcc/system openmpi/4.1.1_ft3 phylobayes-mpi/1.9

cd /mnt/lustre/scratch/nlsas/home/csic/gfy/dri/dzavadska/pb_run/

mpirun -np 64 pb_mpi -s -cat -gtr -d renamedCONCAT.phy cesgaPHYTESTconcat_pb_catgtr_ch2


##############################################
##############################################
##############################################


sbatch -C ilk -t 124:00:00 --mem-per-cpu=1GB cesga_pb_w_renamedCONCAT_ch2.sh

##############################################
##############################################
##############################################
### resuming phylobayes sbatch

#!/bin/bash
#----------------------------------------------------
# PHYLOBAYES SBATCH based on CESGA example 
#----------------------------------------------------
#SBATCH -J cesgaPHYTESTconcat_pb_catgtr_ch1resume            # Job name
#SBATCH -o cesgaPHYTESTconcat_pb_catgtr_ch1res_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e cesgaPHYTESTconcat_pb_catgtr_ch1res_%j.e       # Name of stderr output file(%j expands to jobId)
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB           # Cores per task requested
#SBATCH -t 168:00:00         # Run time (hh:mm:ss) - 2 days     
      
module load cesga/2020 gcc/system openmpi/4.1.1_ft3 phylobayes-mpi/1.9

cd /mnt/lustre/scratch/nlsas/home/csic/gfy/dri/dzavadska/pb_run/

mpirun -np 64 pb_mpi cesgaPHYTESTconcat_pb_catgtr_ch1

##############################################
##############################################
##############################################

sbatch -C ilk -t 168:00:00 --mem-per-cpu=1GB cesga_pb_w_renamedCONCAT_RESUME_ch2.sh
#or
sbatch -t 124:00:00 --mem-per-cpu=1GB cesga_pb_w_renamedCONCAT_RESUME_ch5.sh


###########################################
### tracecomp (if the chains converged)

#on CESGA
module load cesga/2020 gcc/system openmpi/4.1.1_ft3 phylobayes/v4.1e

tracecomp -x 2400 cesgaPHYTESTconcat_pb_catgtr_ch2 cesgaPHYTESTconcat_pb_catgtr_ch4 cesgaPHYTESTconcat_pb_catgtr_ch5 -o ./convergence/pb_run_1200th_gen_Jul16th
#The “-x” option specifies the length of burn-in and the sampling frequency. 
#burn-in of ~20%

bpcomp -x 2400 1 cesgaPHYTESTconcat_pb_catgtr_ch2 cesgaPHYTESTconcat_pb_catgtr_ch4 cesgaPHYTESTconcat_pb_catgtr_ch5 -o ./convergence/pb_run_1200th_gen_Jul16th

bpcomp -x 2400 1 cesgaPHYTESTconcat_pb_catgtr_ch2 cesgaPHYTESTconcat_pb_catgtr_ch4 cesgaPHYTESTconcat_pb_catgtr_ch5 -o ./convergence/pb_run_1200th_gen_Jul9th



#Tracing with Tracer
cd ~/Data/soft/tracer/lib
java -Xmx4096m -jar tracer.jar




#######################################
######### IQTREE SBATCH ###########


#model selection 
iqtree2 -s /scratch5/jpacker/Daryna_proj/actual_trees/CONCAT.fa -m MF -mset WAG,LG,JTT,C10,C20,C30,C40,C50,C60,LG4M,LG4X,CF4 -pre MF_MIX_nopartition_ -safe -nt AUTO

# sbatch ultrafast bootstraps
##############################################
##############################################
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 40

source activate /scratch2/software/anaconda/envs/iqtree-2.2.0.3

iqtree2 -s /scratch5/jpacker/Daryna_proj/actual_trees/CONCAT.fasta -m LG4M+R8 -pre ufboot_LG4Mr8 -bb 1000 -T AUTO

conda deactivate 

##############################################
##############################################

# normal bootstrap 

iqtree2 -s /scratch5/jpacker/Daryna_proj/actual_trees/CONCAT.fasta -m LG4M+R8 -pre realBS250_LG4MR8 -b 250 -T 60

#run stopped on 114th replicate, 
# mapping consensus:
iqtree -t realBS250_LG4MR8.boottrees -sup realBS250_LG4MR8.treefile -pre ML_tree_sup -nt AUTO























################ Code DUMP ############
#Ignore the following


#################################################
# check the overlap with BUSCO dataset

# /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/test_vs_BUSCO
# BUSCO sets downloaded https://busco-data.ezlab.org/v5/data/lineages/




#form a list of all OGs included in the final set and hmms from busco
cd ./euglenozoa_odb10/hmms/
ls * > ../../euglenozoa_odb10_list.txt
cd ../../
cd ./eukaryota_odb10/hmms/
ls * > ../../eukaryota_odb10_list.txt
cd ../../

cat /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/CONCAT.occupancy | cut -f 1 | sed 's/_toMT_gappyout.fasta//g'  > /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/test_vs_BUSCO/all_OGs.txt

ls *_AFTERCHECK.fasta | sed 's/_AFTERCHECK.fasta//g' > /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/test_vs_BUSCO/altered_OGs.txt

#OGs for unaltered

comm -23 <(sort all_OGs.txt) <(sort altered_OGs.txt) > UNaltered_OGs.txt



######## formation of lists finished; the order of the next two lines matters!!
# cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing
# cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/realign 

cat UNaltered_OGs.txt | while read -r og ; do cat euglenozoa_odb10_list.txt | while read -r hmm ; do hmmsearch --tblout ./eugl_outs/"$hmm"_"$og"_eugl ./euglenozoa_odb10/hmms/$hmm /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/$og.fa ; done ; done 

cat altered_OGs.txt | while read -r og ; do cat euglenozoa_odb10_list.txt | while read -r hmm ; do hmmsearch --tblout ./eugl_outs/"$hmm"_"$og"_eugl ./euglenozoa_odb10/hmms/$hmm /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/$og.fa ; done ; done 



cat UNaltered_OGs.txt | while read -r og ; do cat eukaryota_odb10_list.txt | while read -r hmm ; do hmmsearch --tblout ./euk_outs/"$hmm"_"$og"_euk ./eukaryota_odb10/hmms/$hmm /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/$og.fa ; done ; done 


cat altered_OGs.txt | while read -r og ; do cat eukaryota_odb10_list.txt | while read -r hmm ; do hmmsearch --tblout ./euk_outs/"$hmm"_"$og"_euk ./eukaryota_odb10/hmms/$hmm /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/$og.fa ; done ; done 

# master table of all matches

for f in *_eugl; do grep -v "^#" "$f" | sed "s/^/$f\t/" ; done >> all_hits_eugl_outs_odb10.tsv
# "#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------"


for f in *_euk; do grep -v "^#" "$f" | sed "s/^/$f\t/" ; done >> all_hits_euk_outs_odb10.tsv


###### Profile vs profile search


# Build a query HMM from an alignment

cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/*_toMT_gappyout.fasta ./

for f in ./*.fasta ; do hhmake -i $f -M first -o $f.hhm ; done
###
# Combine all .hhm files into a database
# cat *.hhm > ../OG_db.hhm
# Build database index
# hhblitsdb -i OG_db
###


# Create a flatfile + index
ffindex_build -s OG_db.ffdata OG_db.ffindex ./

hhsearch -i ./euglenozoa_odb10/hmms/10092at33682.hmm -d OG_db

ne
for i in $(echo "G65133 RhMon") ; do /home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib "$i" -read1 "$i"*_R1_001.fastq.gz -read2 "$i"*_R2_001.fastq.gz -almosteverything -readlimit 10000000 -readlength 150 ; done




#copy phyloflash outputs

# cp /scratch/data1/agalvez/daryna/phyloflash/G65133.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data1/agalvez/daryna/phyloflash/RhMon.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data2/dzavadska/tr_ass_blastdbs/*.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash

# 

#Download all putative phyloflash contaminant genomes, unzip them and get genomic fnas

cat cont_genomes.txt | while read i ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done 
#for i in $(echo "GCA_030717845.1" "GCA_003259155.1") ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done

#copy all assemblies to denctm folder and create blast databases

cp ./* ../decontam/ass_blastdbs 
cd ../decontam/ass_blastdbs 
for file in *.fasta; do makeblastdb -in $file -dbtype nucl; done


#2) Now BLAST!!!!!!

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query queries.fasta -db /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/ass_blastdbs/"$i"_transcripts.fasta -outfmt 6 -out "$i"_blast_ssu_vs_ass.txt ; done

#manually select putative true matches from the list and create lists of subject contigs to be retrieved and subsequently blasted vs NCBI

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -r -f "$i"_retrieve_cont.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_cont_to_blast.fasta ; done

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query "$i"_cont_to_blast.fasta -db nt -remote -outfmt "6 qseqid sseqid pident length evalue stitle sskingdoms sscinames staxids" -max_target_seqs 1 -out "$i"_blast_contigs_vs_ncbi.txt ; done

#/home/dzavadska/ncbi-blast-2.16.0+/bin/blastn -query D44_cont_to_blast.fasta -db /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -outfmt "6 qseqid sseqid pident length evalue stitle" -max_target_seqs 1 -out D44_blast_contigs_vs_ncbi.txt

#blastall -p blastn -d /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -i test.fasta -o testout.txt


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do perl ./web_blast.pl megablast nt Tabular 1 0 0 "$i"_cont_to_blast.fasta> "$i"_cont_to_blast.fasta ; done
# separately for G65133 because otherwise exceeds memory limit perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_1.fasta > G65133_blast_contigs_vs_ncbi_1.txt
# perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_2.fasta > G65133_blast_contigs_vs_ncbi_2.txt
for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do print "$i" ; grep "Query=" "$i"_blast_contigs_vs_ncbi.txt ; grep "Value" -A 2 "$i"_blast_contigs_vs_ncbi.txt | grep "%" ; done
# grep "Query=" G65133_blast_contigs_vs_ncbi*
# grep "Value" -A 2 G65133_blast_contigs_vs_ncbi* | grep "%"

## manually copy all top hit ids and svae into ncbi_hit_ids.txt

cat ncbi_hit_ids.txt | while read i ; do esummary -db nuccore -id "$i" | xtract -pattern DocumentSummary -element TaxId >> TAXID_ncbi_hit_ids.txt ; done
uniq TAXID_ncbi_hit_ids.txt >> TAXID_ncbi_hit_ids_uniq.txt

#cat TAXID_ncbi_hit_ids_uniq.txt | while read i ; do esearch -db assembly -query ""$i"[TAXID]" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession >> ncbi_blast_cont_genomes.txt ; done
cat TAXID_ncbi_hit_ids_uniq.txt | xargs -n 1 -I{} sh -c 'sleep 1; esearch -db assembly -query "{}[TAXID]" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession' >> ncbi_blast_cont_genomes.txt
mv ncbi_blast_cont_genomes.txt ../contaminant_genomes/

#in ../contaminant_genomes/
cat ncbi_blast_cont_genomes.txt cont_genomes.txt | uniq >> all_cont_genomes_uniq.txt
###Download all contaminant genomes

cat all_cont_genomes_uniq.txt | while read i ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done 

#make blast databases from the downloaded genomes
for i in ./*.fna ; do makeblastdb -in "$i" -dbtype nucl -out "$i" ; done

#BLAST transcriptomes vs contaminant genome databases
/home/drichter/bin/perl/run_BLAST.pl -query /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/ass_blastdbs -database /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/contaminant_genomes -output /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/transcr_vs_cont_genome_out


#subset 60000 random hits
cat ./transcr_vs_cont_genome_out/*decontam_refG_out6.txt > all_decontam_refG.csv
shuf -n 60000 all_decontam_refG.csv > subset_all.csv

#produce a histogram of hits in R, due to huge number of hits, just take 60000 random sample of all hits and plot those using the script in R
# /home/dzavadska/Data/COMPARATIVE_DATASET/phylogen_DECONTAM/blast_res_viz.Rmd


##Remove putative contaminant contigs
for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -v -f "$i"_contaminant_contig_list.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_transcripts_clean.fasta ; done

#double-checking number of contigs before and after decontam

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do echo $i ; grep ">" "$i"_transcripts_clean.fasta | wc -l ; grep ">" ../ass_blastdbs/"$i"_transcripts.fasta | wc -l ; done

# Gemkin 31846 32513
# RhMon 60289 62845
# G65133 209809 248992
# D44 33506 34457
# 16Ckin 33503 33787
# 13G 73868 94536
# 10D 65003 66351

#wcleaner

#wcleaner settings that worked ages ago /scratch/data2/dzavadska/wcleaner_input/settings.yml 

# Batch 1
cp /scratch/data2/dzavadska/wcleaner_input/summer_2023_batch/BEAP* ./

/scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/clean_from_contamination_trass/*_clean.fasta


for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_R1_trimmed.fastq "$i"trimmed_1.fastq ; done
for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_R2_trimmed.fastq "$i"trimmed_2.fastq ; done
for i in $(echo "D44" "16Ckin" "13G" "10D") ; do mv "$i"_transcripts_clean.fasta "$i"trimmed.fasta ; done


for i in $(echo "D44" "16Ckin" "13G" "10D") ; do sed -i 's/_length.*//g' "$i"trimmed.fasta ; done


#Batch 2
fastp --in1 MP2_06939AAF_GCATAGCTAC-AATCGAGGCC_R1_001.fastq.gz --in2 MP2_06939AAF_GCATAGCTAC-AATCGAGGCC_R2_001.fastq.gz --out1 MP2_R1_trimmed.fastq.gz --out2 MP2_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

sudo cp /home/mskamnelou/data/MAOP2_beap0263/MP2_R*_trimmed.fastq.gz

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/fastp/trimmed_D57RNA_06935AAF_AAGAGGATCT-GCTCCTGACG_R*_001.fastq ./
mv trimmed_D57RNA_06935AAF_AAGAGGATCT-GCTCCTGACG_R1_001.fastq D57_R1_trimmed.fastq

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/fastp/trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R*_001.fastq ./
mv trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R1_001.fastq G21_R1_trimmed.fastq
mv trimmed_G21RNA_06936AAF_CCTTCGCGAA-TCGTTCTGCA_R2_001.fastq G21_R2_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/ChY_BEOP0192_Choano_yana_until_trinity/3_Fastp/trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq.gz ./
gunzip trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq.gz
mv trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R1_001.fastq ChY_R1_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/ChY_BEOP0192_Choano_yana_until_trinity/3_Fastp/trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R2_001.fastq ./
mv trimmed_ChY_06937AAF_TTATTCTATG-AACGGACTGC_R2_001.fastq ChY_R2_trimmed.fastq
sudo cp /scratch/data1/mskamnelou/RNAseq/BEAP0192_BEAPO27/BEAP0192_BEAPO279/FRCH_BEAPO279_FRESCHO-3_until_trinity/3_Fastp/trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R*_001.fastq ./
mv trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R1_001.fastq FRCH_R1_trimmed.fastq
mv trimmed_FRCH_06938AAF_GGCACTGCGG-CTCAAGACTT_R2_001.fastq FRCH_R2_trimmed.fastq

sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/Chy_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/d57_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/FRCH_1st_decont_trinity_out_dir.Trinity.fasta ./
sudo cp /scratch/data1/agalvez/transcriptomes/2nd_batch/trinity/1st_decont/output/clean/g21_1st_decont_trinity_out_dir.Trinity.fasta ./


for i in $(echo "RhMon" "G65133" "MP2") ; do sed -i 's/_length.*//g' "$i"trimmed.fasta ; done
for i in $(echo "G21" "ChY" "FRCH" "D57") ; do sed -i 's/ len.*//g' "$i"trimmed.fasta ; done


###############################
# CROSS-CONTAMINATION REMOVAL #
###############################
#RUNNNING WCLEANER

#Run data prep; takes some time ~o/n
winston_prepare_data
#Run cleaner; takes 0 time
winston_find_contaminations

#deactivate wcleaner 
deactivate

#WCLEANER OUT

#Check the number of deleted contigs from *deleted.fasta files

for file in *deleted.fasta ; do echo $file ; grep '>' $file | wc -l ; done

# 10Dtrimmed_deleted.fasta
# 325
# 13Gtrimmed_deleted.fasta
# 915
# 16Cinittrimmed_deleted.fasta
# 56728
# 16Ckintrimmed_deleted.fasta
# 712
# D44trimmed_deleted.fasta
# 166

# G65133trimmed_deleted.fasta
# 146
# RhMontrimmed_deleted.fasta
# 2771


#check distribution of contig length (after decontam) - just in case
awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' Bedax.fasta > Bedax_len_distr.txt
awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' Neocurv.fasta > Neocurv_len_distr.txt

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' "$i"_doubleclean_contigs.fasta > "$i"_len_distr.txt ; done





#########################
# PROTEIN PREDICTION 
#########################



###########################
#Running TransDecoder


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do mv "$i"*.fasta "$i"_doubleclean_contigs.fasta ; done
cp Gemkin_doubleclean_contigs.fasta Gemkin_doubleclean_contigs_RC.fasta

for file in ./Gemkin*_doubleclean_contigs_RC.fasta ; do sudo TransDecoder.LongOrfs -t $file ; done
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file - therefor, ThransDecoder ran in a different mode

for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.LongOrfs -t $file -S ; done


for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.Predict -t $file ; done

#########################
# BUSCO COMPLETENESS 
#########################

#running BUSCO

# cp ./*transdecoder.pep ./busco
# cp ./*transdecoder.cds ./busco

for file in *transdecoder.pep ; do busco -i $file -l eukaryota_odb10 -o "$file"pep -m proteins --cpu 8 ; done
for file in *transdecoder.pep ; do busco -i $file -l euglenozoa_odb10 -o "$file"EUGLENOZOApep -m proteins --cpu 8 ; done
















####################################################################
####################################################################
####################################################################
# Orthology inference and verification
####################################################################
####################################################################
####################################################################


# OrthoFinder 45 taxa, CD-hit 90%, diamond_ultra_sens , msa

/home/agalvez/programs/OrthoFinder/orthofinder -f NR_proteins_90_45taxa -S diamond_ultra_sens -M msa -n "Results_diamond_ultra_sens_FASTTREE" -a 16 -t 16


#to obtain statistics per species 

# cp ./SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv  ./Summaries_SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa_Results_diamond_ultra_sens_FASTTREE_Statistics_PerSpecies.tsv

# cp ./SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Orthogroups/Orthogroups.GeneCount.tsv ./Summaries_SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa_Results_diamond_ultra_sens_FASTTREE_Orthogroups.GeneCount.tsv




## to produce "OG_summary_NR_proteins_45_taxa_diamond_ultra_sens.csv" , run Orthofinder_stats_UNIVERSAL.Rmd

#filtering in R from OG_summary:

#OG_local <- fread("OG_summary_NR_proteins_45_taxa_diamond_ultra_sens.csv")
#subset <- OG_local[which(OG_local$mean<=1 & OG_local$sd<=1 & OG_local$sp_repr>25),]
#double_subset <- subset[which(subset$mean>=0.9 & subset$sd<=0.5 & subset$sp_repr>35),]
#fwrite(double_subset, "double_subset_NR_proteins_45_taxa_diamond_ultra_sens.csv")

## 494 OG names saved into OG_double_subset.txt


#cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing
cat OG_double_subset.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Gene_Trees/"$i"* ./ ; done

cat OG_double_subset.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/NR_proteins_90_45taxa/OrthoFinder/Results_diamond_ultra_sens_FASTTREE/Orthogroup_Sequences/"$i"* ./ ; done

# renaming sequence headers

cat OG_double_subset.txt | while read i ; do sed 's/;/_/g' ./"$i".fa > ./"$i"_charsubs.fa ; done
cat OG_double_subset.txt | while read i ; do phykit tip_labels ./"$i"_tree.txt > ./"$i"_tiplabels.txt ; done
cat OG_double_subset.txt | while read i ; do grep ">" ./"$i"_charsubs.fa | sed 's/>//g' > ./"$i"_fastalabels.txt ; done
cat OG_double_subset.txt | while read e ; do while read i; do sed -i "s@"$i"@$(grep -F "$i" ./"$e"_tiplabels.txt)@g" ./"$e"_charsubs.fa ; done < "$e"_fastalabels.txt ; done
cat OG_double_subset.txt | while read e ; do mv ./"$e"_charsubs.fa ./"$e"_sp_renamed.fasta ; done


#aligning, trimming and building tree for SG-tree check
cat OG_double_subset.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 16 --reorder ""$i"_sp_renamed.fasta"  > ""$i"_sp_renamed_MAFFT.fasta" ; done

#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/TRIMMED_ASSESSED
cat ./OG_double_subset.txt | while read i ; do trimal -in "$i"_sp_renamed_MAFFT.fasta -out "$i"_mafft_gappyout.fasta -gappyout -fasta ; done

#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/IQTREE_sg/gappyout_iqtree/
for i in ./*_mafft_gappyout.fasta ; do iqtree -bb 1000 -seed 1234 -nt 24 -s "$i" ; done
# generates the set of trees used FOR THE ACTUAL CHECKS; *.contree


##### After the manual check

#copy trees and alignments that were checked from dazhbog  
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/downstream_orthogroups/orthosnap/Genes_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/downstream_orthogroups/orthosnap/Trees_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN
#copy the initial OG sequences to the new directory
# cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/*_sp_renamed.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/TO_CHECK_UPDATED_no_orthoSNAP/

#replacing | in fasta headers with _
for i in ./TO_CHECK_UPDATED_no_orthoSNAP/* ; do sed -i 's/|/_/g' $i ; done

#sed 's/|/_/g' ./removing_genes/OG0003302_sp_renamed.fasta #debug, check Sspe gene names
#check 
grep -F "|" ./TO_CHECK_UPDATED_no_orthoSNAP/*


# deleting genes

cat Genes_to_delete.tsv | cut -f 1 | uniq > uniq_OG_list.txt

cat uniq_OG_list.txt | while read i ; do grep "$i" Genes_to_delete.tsv | cut -f 2 > ./removing_genes/"$i"_to_delete.txt ; done

cat uniq_OG_list.txt | while read i ; do seqkit grep -n -v -f ./removing_genes/"$i"_to_delete.txt ./TO_CHECK_UPDATED_no_orthoSNAP/"$i"_sp_renamed.fasta -o "$i"_AFTERCHECK.fasta ; done


## Check if everything was removed
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep "$e" "$i"_AFTERCHECK.fasta ; echo "$i" ; done ; done
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep -F "$e" "$i"_AFTERCHECK.fasta ; done ; done

#manually remove few sequences that failed to be removed by seqkit for whatever reason
#>RhMon_Gene.33090__NODE_18304__g.33090__m.33090
#OG0002451
#>Cruzella_marina_Gene.84135__NODE_19292_length_2400_cov_477.822643_g7898_i0__g.84135__m.84135
#OG0002445
#>G65133_Gene.62375__NODE_45689__g.62375__m.62375
#OG0002452
#>RhMon_Gene.25779__NODE_13194__g.25779__m.25779
#OG0002452


# seqkit grep -v -f ./removing_genes/OG0003302_to_delete.txt ./TO_CHECK_UPDATED_no_orthoSNAP/OG0003302_sp_renamed.fasta

########## REALIGNING

mv ./*_AFTERCHECK.fasta ./realign/
cd ./realign/
cat ../uniq_OG_list.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 100 --reorder ""$i"_AFTERCHECK.fasta"  > ""$i"_AFTERCHECK_realigned.fasta" ; done



cat ./uniq_OG_list.txt | while read i ; do cp ./realign/"$i"_AFTERCHECK_realigned.fasta ./to_modeltest/ ; done

## removing trees marked for complete removal
cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest
cat ../Trees_to_delete.tsv | while read i ; do rm ./"$i"_AFTERCHECK_realigned.fasta ; done


#trimming the re-aligned ones

#for i in ./*AFTERCHECK_realigned.fasta  ; do trimal -in "$i" -out "$i"_gappyout.fasta -gappyout -fasta ; done
#WARNING: Removing sequence 'Gene.47917__CAMNT_0016070987__g.47917__m.47917' composed only by gaps
#WARNING: Removing sequence 'Gene.34531__NODE_19382__g.34531__m.34531' composed only by gaps

cat ./Genes_to_delete.tsv | cut -f 1 | sort -u | > OGs_realigned.txt

##adding alignments which are meant to be kept as they were initially
cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/*_sp_renamed_MAFFT.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/tmp/
cat ./Trees_to_delete.tsv | while read i ; do rm ./tmp/"$i"* ; done
cat ./uniq_OG_list.txt | while read i ; do rm ./tmp/"$i"* ; done
cp ./tmp/* ./to_modeltest/



#
ls | sed 's/_.*//g' > ../All_OGs_kept.txt
#
cat ../All_OGs_kept.txt | while read i ; do trimal -in ./"$i"* -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done

#WARNING: Removing sequence 'G65133_Gene.88879__NODE_82234__g.88879__m.88879' composed only by gaps

#renaming for concatenation

#cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/FASTTREE_draft/gappyout_fasttree/species_names.txt /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/
for file in ./*_gappyout.fasta ; do cat ../species_names.txt | while read i ; do sed -i "s/${i}_/${i}|/g" "$file" ; done ; done

## checking if there are duplicates
for i in ./*_gappyout.fasta ; do echo $(grep ">" "$i" | sed 's/|.*//g' | sort | uniq -cd ) ; echo "$i" ; done

## checking if every header got substituted with | species separator
for i in ./*_gappyout.fasta ; do echo "$i $(grep -c '^>' "$i") $(grep -c '|' "$i")" ; done

#re-doing a few duplicates detected, just manually on dazhbog
#2 >G65133
#./OG0002646_toMT_gappyout.fasta
#2 >Strigomonas_culicis
#./OG0002823_toMT_gappyout.fasta
#2 >G65133
#./OG0002831_toMT_gappyout.fasta

scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/OG0002831* ./
for i in ./* ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 100 --reorder "$i"  > ""$i"_AFTERCHECK_realigned.fasta" ; done
for i in $(echo "OG0002646" "OG0002823" "OG0002831")  ; do /home/dzavadska/Data/soft/trimAl/source/trimal -in ./"$i"_AFTERCHECK_realigned.fasta_AFTERCHECK_realigned.fasta -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done
scp ./*gappyout.fasta nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest/



# for i in $(echo "OG0002650" "OG0003130" "OG0003260" "OG0003298")  ; do trimal -in ./"$i"*_realigned.fasta -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done
for file in $(echo  "OG0002646" "OG0002823" "OG0002831") ; do cat ../species_names.txt | while read i ; do sed -i "s/${i}_/${i}|/g" "$file"*_gappyout.fasta ; done ; done

#substitute specifically >g31 to Willaertia_magna|g31 - for some reason I dont really remeber, this is some strange exception.
grep ">g31" ./*_gappyout.fasta
#./OG0002815_toMT_gappyout.fasta:>g31
sed -i 's/>g31/>Willaertia_magna|g31/g' ./OG0002815_toMT_gappyout.fasta
grep ">" ./OG0002815_toMT_gappyout.fasta

#####################################################################################################
############################# CONCATENATING AND TESTING MODELS ######################################
#####################################################################################################
ls *_gappyout.fasta > alignments.txt

#replacing | with " "
for i in *_gappyout.fasta ; do cp "$i" ./ "$i"_BACKUP.fasta ; done
# debug
#sed 's/|/ /g' OG0003327_toMT_gappyout.fasta_BACKUP.fasta | head
for i in *_gappyout.fasta ; do sed -i 's/|/ /g' "$i" ; done

pk_create_concat -a alignments.txt -p CONCAT

#check if concatenation was done correctly; there should be only the names of taxa, same number as that of intital taxa in a dataset
grep ">" CONCAT.fa

#Partition file output: CONCAT.partition
#Concatenated fasta output: CONCAT.fa
#Occupancy report: CONCAT.occupancy

cp ./CONCAT* ../../../CONCAT_AND_GIANT_TREES/model_selection/








