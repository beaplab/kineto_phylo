
#########################
# FASTQC (sanity check)
#########################
fastqc --noextract --nogroup -o fastqc *.fastq.gz

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
AGATCGGAA



###batch 2025
for i in $(echo "BABKIN KUCKIN INUKIN LUKPR6 Gre177 MP2400 Pmic2 PN") ; do fastp --in1 "$i"*_1.fastq.gz --in2 "$i"*_2.fastq.gz --out1 "$i"_R1_trimmed.fastq.gz --out2 "$i"_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ; done


fastp --in1 OB23*_1.fastq.gz --in2 OB23*_2.fastq.gz --out1 OB23_R1_trimmed.fastq.gz --out2 OB23_R2_trimmed.fastq.gz --low_complexity_filter --cut_front --cut_tail --cut_right --cut_front_window_size 1 --cut_tail_window_size 1 --cut_right_window_size 4 --cut_mean_quality 5 --trim_front1 12 --trim_front2 12 --trim_poly_g --trim_poly_x --adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG


#########################
# ASSEMBLY 
#########################
# (originating from different batches, reads from different cultures were assembles at different times, although with technically similar commands. For the record, I keep all of them as they were ran, and the preceeding line corresponds to the name of the script which contains original commands ran on a given batch)

#trassembly_phylogenomics.sh
# in /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/ on nisaba
rnaspades -1 GK1_R1_trimmed.fastq.gz -2 GK1_R2_trimmed.fastq.gz --ss-fr -o /scratch/data2/dzavadska/RNAspades/Gemkin/
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file 


#alt_assemblers_trassembly_auto_attempt.sh - original source of the following line
for i in $(echo "D44 10D 13G 16Ckin") ; do nohup rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done


#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do rnaspades -1 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/transcriptome_reads/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/
 ; done


###batch 2025
for i in $(echo "BABKIN KUCKIN INUKIN LUKPR6 Gre177 MP2400 Pmic2 PN") ; do rnaspades -1 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/"$i"_R2_trimmed.fastq.gz --ss-rf -o /scratch/data2/dzavadska/RNAspades/"$i"/ ; done

# in /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/ on nisaba
rnaspades -1 OB23_R1_trimmed.fastq.gz -2 OB23_R2_trimmed.fastq.gz -o /scratch/data2/dzavadska/RNAspades/OB23/





#########################
# DECONTAMINATION 
#########################

###!!!!!!!!!!!! for 2025 batch, similar procedures were replicated in /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/, filenames and loops work universally!


##############
# PhyloFlash #
##############

#trassembly_phylogenomics.sh
/home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib Gemkin -read1 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R1_trimmed.fastq.gz -read2 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/GK1_R2_trimmed.fastq.gz -almosteverything -readlimit 10000000 -readlength 86

## for "D44 10D 13G 16Ckin" , PhyloFlash was ran by Daniel

#trassembly_Jun2024_auto_attempt.sh - original source of the following line
for i in $(echo "G65133 RhMon") ; do /home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib "$i" -read1 "$i"*_R1_001.fastq.gz -read2 "$i"*_R2_001.fastq.gz -almosteverything -readlimit 10000000 -readlength 150 ; done

## for "BABKIN KUCKIN INUKIN LUKPR6 Gre177 MP2400 Pmic2 PN" , PhyloFlash was ran by Daniel

/home/dzavadska/lib/phyloFlash/phyloFlash.pl -lib "OB23" -read1 OB23*_1.fastq.gz -read2 OB23*_2.fastq.gz -almosteverything -readlimit 10000000 -readlength 150



#copy phyloflash outputs

# cp /scratch/data1/agalvez/daryna/phyloflash/G65133.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data1/agalvez/daryna/phyloflash/RhMon.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash
# cp /scratch/data2/dzavadska/tr_ass_blastdbs/*.spades_rRNAs.final.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/phyloflash


# !!! ALL transcriptome PhyloFlash htmls on dazhbog in /home/dzavadska/Data/transcriptomes/phyloFlash/


# save contaminant genome IDs to cont_genomes.txt

#Download all putative phyloflash contaminant genomes, unzip them and get genomic fnas

cat cont_genomes.txt | while read i ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done 
#for i in $(echo "GCA_030717845.1" "GCA_003259155.1") ; do datasets download genome accession $i --include genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/* ./ ; rm -r ./ncbi_dataset ; done

#copy all assemblies to denctm folder and create blast databases

cp ./* ../decontam/ass_blastdbs 
cd ../decontam/ass_blastdbs 
for file in *.fasta; do makeblastdb -in $file -dbtype nucl; done

# batch 2
# for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do cp /scratch/data2/dzavadska/RNAspades/"$i"/transcripts.fasta ./ass_blastdbs/"$i"_transcripts.fasta ; done 


#2) Now BLAST!!!!!!

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query queries.fasta -db /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam/ass_blastdbs/"$i"_transcripts.fasta -outfmt 6 -out "$i"_blast_ssu_vs_ass.txt ; done

for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do blastn -query queries.fasta -db /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/ass_blastdbs/"$i"_transcripts.fasta -outfmt 6 -out "$i"_blast_ssu_vs_ass.txt ; done

#manually select putative true matches from the list and create lists of subject contigs to be retrieved and subsequently blasted vs NCBI

# save the list of putative contaminant 18S contigs to *_retrieve_cont.txt


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -r -f "$i"_retrieve_cont.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_cont_to_blast.fasta ; done

for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do seqkit grep -r -f "$i"_retrieve_cont.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_cont_to_blast.fasta ; done

# the commented versions below hit ncbi timeout; used web_blast.pl by Daniel, which implements megablast
#for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do blastn -query "$i"_cont_to_blast.fasta -db nt -remote -outfmt "6 qseqid sseqid pident length evalue stitle sskingdoms sscinames staxids" -max_target_seqs 1 -out "$i"_blast_contigs_vs_ncbi.txt ; done

#/home/dzavadska/ncbi-blast-2.16.0+/bin/blastn -query D44_cont_to_blast.fasta -db /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -outfmt "6 qseqid sseqid pident length evalue stitle" -max_target_seqs 1 -out D44_blast_contigs_vs_ncbi.txt

#blastall -p blastn -d /scratch/datasets_symbolic_links/common_databases/nt/nt_blastdb/nt_blast -i test.fasta -o testout.txt


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do perl ./web_blast.pl megablast nt Tabular 1 0 0 "$i"_cont_to_blast.fasta > "$i"_blast_contigs_vs_ncbi.txt  ; done
# separately for G65133 because otherwise exceeds memory limit perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_1.fasta > G65133_blast_contigs_vs_ncbi_1.txt
# perl ./web_blast.pl megablast nt Tabular 1 0 0 G65133_cont_to_blast_2.fasta > G65133_blast_contigs_vs_ncbi_2.txt

for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do perl ./web_blast.pl megablast nt Tabular 1 0 0 "$i"_cont_to_blast.fasta > "$i"_blast_contigs_vs_ncbi.txt  ; done

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do print "$i" ; grep "Query=" "$i"_blast_contigs_vs_ncbi.txt ; grep "Value" -A 2 "$i"_blast_contigs_vs_ncbi.txt | grep "%" ; done
# grep "Query=" G65133_blast_contigs_vs_ncbi*
# grep "Value" -A 2 G65133_blast_contigs_vs_ncbi* | grep "%"
for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do print "$i" ; grep "Query=" "$i"_blast_contigs_vs_ncbi.txt ; grep "Value" -A 2 "$i"_blast_contigs_vs_ncbi.txt | grep "%" ; done

## manually copy all top hit ids and svae into ncbi_hit_ids.txt


cat ncbi_hit_ids.txt | while read i ; do esummary -db nuccore -id "$i" | xtract -pattern DocumentSummary -element TaxId >> TAXID_ncbi_hit_ids.txt ; done
sort -u TAXID_ncbi_hit_ids.txt > TAXID_ncbi_hit_ids_uniq.txt



# removing all unclutured organisms from the list  
#form the list of taxid AND Organism name
# cat ncbi_hit_ids.txt | while read i ; do esummary -db nuccore -id "$i" | xtract -pattern DocumentSummary -element TaxId Organism >> Organism_ncbi_hit_ids.txt ; done
# uniq Organism_ncbi_hit_ids.txt >> Organism_ncbi_hit_ids_uniq.txt
# grep "uncultured" Organism_ncbi_hit_ids_uniq.txt | cut -f 1 | sort -u > uncultured_cont_txids.txt
# manually add E.coli taid to the list : 562 !!!
# cat uncultured_cont_txids.txt | while read i ; do sed -i "s/$i//g" TAXID_ncbi_hit_ids_uniq.txt ; done

# grep "77133" TAXID_ncbi_hit_ids_uniq.txt # sanity check
# grep "562" TAXID_ncbi_hit_ids_uniq.txt # sanity check



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


#BLAST transcriptomes vs contaminant genome databases
/home/drichter/bin/perl/run_BLAST.pl -query /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/ass_blastdbs -database /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/contaminant_genomes -output /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/transcr_vs_cont_genome_out



#subset 60000 random hits
cat ./transcr_vs_cont_genome_out/*decontam_refG_out6.txt > all_decontam_refG.csv
shuf -n 60000 all_decontam_refG.csv > subset_all.csv

#produce a histogram of hits in R, due to huge number of hits, just take 60000 random sample of all hits and plot those using the script in R
# /home/dzavadska/Data/COMPARATIVE_DATASET/phylogen_DECONTAM/blast_res_viz.Rmd


_contaminant_contig_list.txt


##Remove putative contaminant contigs
for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do seqkit grep -v -f "$i"_contaminant_contig_list.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_transcripts_clean.fasta ; done

for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do seqkit grep -v -f "$i"_contaminant_contig_list.txt ../ass_blastdbs/"$i"_transcripts.fasta >> "$i"_transcripts_clean.fasta ; done



#double-checking number of contigs before and after decontam

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do echo $i ; grep ">" "$i"_transcripts_clean.fasta | wc -l ; grep ">" ../ass_blastdbs/"$i"_transcripts.fasta | wc -l ; done

# Gemkin 31846 32513
# RhMon 60289 62845
# G65133 209809 248992
# D44 33506 34457
# 16Ckin 33503 33787
# 13G 73868 94536
# 10D 65003 66351


for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do echo $i ; grep ">" "$i"_transcripts_clean.fasta | wc -l ; grep ">" ../ass_blastdbs/"$i"_transcripts.fasta | wc -l ; done

# OB23 42526 43440
# BABKIN 40074 41319
# KUCKIN 85700 86886
# INUKIN 27886 28122
# LUKPR6 34804 39557
# Gre177 21189 21248




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



#batch 2.2 (aug2025)
#in /scratch/data2/dzavadska/COMPARATIVE_DATASET/decontam_batch2/wcleaner/batch2_2_input

cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/*trimmed* ./

mv BABKIN_R1_trimmed.fastq.gz BABKINtrimmed_1.fastq.gz
mv BABKIN_R2_trimmed.fastq.gz BABKINtrimmed_2.fastq.gz
mv INUKIN_R1_trimmed.fastq.gz INUKINtrimmed_1.fastq.gz
mv INUKIN_R2_trimmed.fastq.gz INUKINtrimmed_2.fastq.gz
mv KUCKIN_R1_trimmed.fastq.gz KUCKINtrimmed_1.fastq.gz
mv KUCKIN_R2_trimmed.fastq.gz KUCKINtrimmed_2.fastq.gz
mv LUKPR6_R1_trimmed.fastq.gz LUKPR6trimmed_1.fastq.gz
mv LUKPR6_R2_trimmed.fastq.gz LUKPR6trimmed_2.fastq.gz
mv OB23_R1_trimmed.fastq.gz OB23trimmed_1.fastq.gz
mv OB23_R2_trimmed.fastq.gz OB23trimmed_2.fastq.gz
mv Gre177_R1_trimmed.fastq.gz Gre177trimmed_1.fastq.gz
mv Gre177_R2_trimmed.fastq.gz Gre177trimmed_2.fastq.gz

mv MP2400_R1_trimmed.fastq.gz MP2400trimmed_1.fastq.gz
mv MP2400_R2_trimmed.fastq.gz MP2400trimmed_2.fastq.gz
mv PN_R1_trimmed.fastq.gz PNtrimmed_1.fastq.gz
mv PN_R2_trimmed.fastq.gz PNtrimmed_2.fastq.gz
mv Pmic2_R1_trimmed.fastq.gz Pmic2trimmed_1.fastq.gz
mv Pmic2_R2_trimmed.fastq.gz Pmic2trimmed_2.fastq.gz


for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177 MP2400 Pmic2 PN") ; do cp /scratch/data2/dzavadska/RNAspades/"$i"/transcripts.fasta ./"$i"trimmed.fasta ; done 

# VERY IMPORTANT !!!!!!!
for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177 MP2400 Pmic2 PN") ; do sed -i 's/_length.*//g' "$i"trimmed.fasta ; done
#sanity check afterwards
# grep ">" LUKPR6trimmed.fasta | wc -l
# grep ">" LUKPR6trimmed.fasta | sort -u | wc -l
 



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


# BABKINtrimmed_deleted.fasta
# 39
# Gre177trimmed_deleted.fasta
# 49
# INUKINtrimmed_deleted.fasta
# 95
# KUCKINtrimmed_deleted.fasta
# 6
# LUKPR6trimmed_deleted.fasta
# 25471
# MP2400trimmed_deleted.fasta
# 12
# Pmic2trimmed_deleted.fasta
# 101
# PNtrimmed_deleted.fasta
# 65


#check distribution of contig length (after decontam) - just in case

for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' "$i"_doubleclean_contigs.fasta > "$i"_len_distr.txt ; done



#########################
# PROTEIN PREDICTION 
#########################



###########################
#Running TransDecoder


for i in $(echo "Gemkin" "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do mv "$i"*.fasta "$i"_doubleclean_contigs.fasta ; done
for i in $(echo "OB23 BABKIN KUCKIN INUKIN LUKPR6 Gre177") ; do mv "$i"*.fasta "$i"_doubleclean_contigs.fasta ; done
cp Gemkin_doubleclean_contigs.fasta Gemkin_doubleclean_contigs_RC.fasta



for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.LongOrfs -t $file -S ; done
for file in ./*_doubleclean_contigs.fasta ; do sudo TransDecoder.Predict -t $file ; done



for file in ./Gemkin*_doubleclean_contigs_RC.fasta ; do sudo TransDecoder.LongOrfs -t $file ; done
#*** here the FR orientation was in fact incorrectly specified - look up the end of the file - therefor, ThransDecoder ran in a different mode
for file in ./OB23*_doubleclean_contigs.fasta ; do sudo TransDecoder.LongOrfs -t $file ; done
for file in ./OB23*_doubleclean_contigs.fasta ; do sudo TransDecoder.Predict -t $file ; done
#non-strand-specific sequencing platform

#########################
# BUSCO COMPLETENESS 
#########################

#running BUSCO

# cp ./*transdecoder.pep ./busco
# cp ./*transdecoder.cds ./busco

for file in *transdecoder.pep ; do busco -i $file -l eukaryota_odb10 -o "$file"pep -m proteins --cpu 8 ; done
for file in *transdecoder.pep ; do busco -i $file -l euglenozoa_odb10 -o "$file"EUGLENOZOApep -m proteins --cpu 8 ; done



#############################
# CD-HIT REDUNDANCY REDUCTION
#############################

for i in $(echo "OB23 BABKIN KUCKIN INUKIN Gre177") ; do cd-hit -i ./"$i"*transdecoder.pep -o "$i"_NR90_prot.fasta ; done











####################################################################
####################################################################
####################################################################
# Orthology inference and verification
####################################################################
####################################################################
####################################################################


# OrthoFinder 49 taxa, CD-hit 90%, diamond_ultra_sens , msa

#here /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa
/home/agalvez/programs/OrthoFinder/orthofinder -f NR_proteins_90_49taxa -S diamond_ultra_sens -M msa -n "Results_diamond_ultra_sens_FASTTREE1" -a 16 -t 16



#on dazhbog

# cd /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY

# scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa/OrthoFinder/Results_Results_diamond_ultra_sens_FASTTREE1/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv ./

# scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa/OrthoFinder/Results_Results_diamond_ultra_sens_FASTTREE1/Orthogroups/Orthogroups.GeneCount.tsv ./




## to produce "OG_summary_V2.csv" , run Orthofinder_stats_UNIVERSAL_V2.Rmd

#filtering in R from OG_summary:

# OG_local <- fread("OG_summary_V2.csv")
# subset <- OG_local[which(OG_local$mean<=1 & OG_local$sd<=1 & OG_local$sp_repr>25),]
#double_subset <- subset[which(subset$mean>=0.9 & subset$sd<=0.5 & subset$sp_repr>35),]
#fwrite(double_subset, "double_subset_V2.csv")

## 461 OG names saved into OG_double_subsetV2.txt

#cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/orthogroup_postprocessing
cat OG_double_subsetV2.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa/OrthoFinder/Results_Results_diamond_ultra_sens_FASTTREE1/Gene_Trees/"$i"* ./ ; done

cat OG_double_subsetV2.txt | while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa/OrthoFinder/Results_Results_diamond_ultra_sens_FASTTREE1/Orthogroup_Sequences/"$i"* ./ ; done


# renaming sequence headers

cat OG_double_subsetV2.txt | while read i ; do sed 's/;/_/g' ./"$i".fa > ./"$i"_charsubs.fa ; done
cat OG_double_subsetV2.txt | while read i ; do phykit tip_labels ./"$i"_tree.txt > ./"$i"_tiplabels.txt ; done

cat OG_double_subsetV2.txt | while read i ; do grep ">" ./"$i"_charsubs.fa | sed 's/>//g' > ./"$i"_fastalabels.txt ; done
cat OG_double_subsetV2.txt | while read e ; do while read i; do sed -i "s@"$i"@$(grep -F "$i" ./"$e"_tiplabels.txt)@g" ./"$e"_charsubs.fa ; done < "$e"_fastalabels.txt ; done
cat OG_double_subsetV2.txt | while read e ; do mv ./"$e"_charsubs.fa ./"$e"_sp_renamed.fasta ; done




#aligning, trimming and building tree for SG-tree check
cat OG_double_subsetV2.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 16 --reorder ""$i"_sp_renamed.fasta"  > ""$i"_sp_renamed_MAFFT.fasta" ; done

#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/TRIMMED_ASSESSED
cat ../OG_double_subsetV2.txt | while read i ; do trimal -in "$i"_sp_renamed_MAFFT.fasta -out "$i"_mafft_gappyout.fasta -gappyout -fasta ; done
# WARNING: Removing sequence 'RhMon_Gene.34531__NODE_19382__g.34531__m.34531' composed only by gaps
# WARNING: Removing sequence 'Ophirina_chinija_NODE_19866_length_524_cov_5.328859_g11489_i0.p2' composed only by gaps
# WARNING: Removing sequence 'D44_Gene.55280__NODE_26646__g.55280__m.55280' composed only by gaps
# WARNING: Removing sequence 'OB23_Gene.104073__NODE_25617__g.104073__m.104073' composed only by gaps



#/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/IQTREE_sg/gappyout_iqtree/
for i in ./*_mafft_gappyout.fasta ; do iqtree -bb 1000 -seed 1234 -nt 24 -s "$i" ; done
# generates the set of trees used FOR THE ACTUAL CHECKS; *.contree





##### After the manual check

#copy trees and alignments that were checked from dazhbog  
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/trees_to_check/Genes_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/MANUAL_TREE_CURATION_REALIGN
# scp /home/dzavadska/Data/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/trees_to_check/Trees_to_delete.tsv  nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/MANUAL_TREE_CURATION_REALIGN
#copy the initial OG sequences to the new directory
# cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/orthogroup_postprocessing/*_sp_renamed.fasta /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/MANUAL_TREE_CURATION_REALIGN/seqs

#replacing | in fasta headers with _
for i in ./seqs/* ; do sed -i 's/|/_/g' $i ; done

#sed 's/|/_/g' ./removing_genes/OG0003302_sp_renamed.fasta #debug, check Sspe gene names
#check 
grep -F "|" ./seqs/*


# deleting genes

cat Genes_to_delete.tsv | cut -f 1 | uniq > uniq_OG_list.txt

# mkdir removing_genes
cat uniq_OG_list.txt | while read i ; do grep "$i" Genes_to_delete.tsv | cut -f 2 > ./removing_genes/"$i"_to_delete.txt ; done

cat uniq_OG_list.txt | while read i ; do seqkit grep -n -v -f ./removing_genes/"$i"_to_delete.txt ./seqs/"$i"_sp_renamed.fasta -o "$i"_AFTERCHECK.fasta ; done


## Check if everything was removed
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep "$e" "$i"_AFTERCHECK.fasta ; echo "$i" ; done ; done
cat uniq_OG_list.txt | while read i ; do cat ./removing_genes/"$i"_to_delete.txt | while read e ; do grep -F "$e" "$i"_AFTERCHECK.fasta ; done ; done


#manually remove few sequences that failed to be removed by seqkit for whatever reason - in this case, none detected

#moving OGs left intact 
cat ./tokeep_OG_list.txt |  while read i ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/orthogroup_postprocessing/"$i"_sp_renamed.fasta ./realign/"$i"_AFTERCHECK.fasta ; done




#cat ./uniq_OG_list.txt |  while read i ; do grep "$i" tokeep_OG_list.txt ; done
#cat all_OGs.txt |  while read i ; do ls ./realign/"$i"_AFTERCHECK.fasta ; done

########## REALIGNING
#mkdir realign

mv ./*_AFTERCHECK.fasta ./realign/
cd ./realign/
cat ../uniq_OG_list.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 1000 --reorder ""$i"_AFTERCHECK.fasta"  > ""$i"_AFTERCHECK_realigned.fasta" ; done
cat ../tokeep_OG_list.txt | while read i ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 1000 --reorder ""$i"_AFTERCHECK.fasta"  > ""$i"_AFTERCHECK_realigned.fasta" ; done
# for i in $(echo "OG0002773" "OG0003014") ; do "/usr/bin/mafft" --thread 24 --localpair --maxiterate 1000 --reorder ""$i"_AFTERCHECK.fasta"  > ""$i"_AFTERCHECK_realigned.fasta" ; done

cat ./uniq_OG_list.txt | while read i ; do cp ./realign/"$i"_AFTERCHECK_realigned.fasta ./to_modeltest/ ; done
cat ./tokeep_OG_list.txt | while read i ; do cp ./realign/"$i"_AFTERCHECK_realigned.fasta ./to_modeltest/ ; done
# for i in $(echo "OG0002773" "OG0003014") ; do cp ./realign/"$i"_AFTERCHECK_realigned.fasta ./to_modeltest/ ; done


## removing trees marked for complete removal (they must not be there at all, to begin with, but just in case)
# cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/MANUAL_TREE_CURATION_REALIGN/to_modeltest
# cat ../Trees_to_delete.tsv | while read i ; do rm ./"$i"_AFTERCHECK_realigned.fasta ; done





#forming list of all OGs kept
ls | sed 's/_.*//g' > ../All_OGs_kept.txt

#trimming the re-aligned ones
cat ../All_OGs_kept.txt | while read i ; do trimal -in ./"$i"_AFTERCHECK_realigned.fasta -out "$i"_toMT_gappyout.fasta -gappyout -fasta ; done

# WARNING: Removing sequence 'RhMon_Gene.34531__NODE_19382__g.34531__m.34531' composed only by gaps
# WARNING: Removing sequence 'OB23_Gene.104073__NODE_25617__g.104073__m.104073' composed only by gaps
# WARNING: Removing sequence 'D44_Gene.55280__NODE_26646__g.55280__m.55280' composed only by gaps


#renaming for concatenation

#cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_1ST_TRY/orthogroup_postprocessing/FASTTREE_draft/gappyout_fasttree/species_names.txt /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/MANUAL_TREE_CURATION_REALIGN/
#MANUALLY ADD OB23, KUCKIN, INUKIN, BABKIN
for file in ./*_gappyout.fasta ; do cat ../species_names.txt | while read i ; do sed -i "s/${i}_/${i}|/g" "$file" ; done ; done

## checking if there are duplicates
for i in ./*_gappyout.fasta ; do echo $(grep ">" "$i" | sed 's/|.*//g' | sort | uniq -cd ) ; echo "$i" ; done

## checking if every header got substituted with | species separator
for i in ./*_gappyout.fasta ; do echo "$i $(grep -c '^>' "$i") $(grep -c '|' "$i")" ; done


#substitute specifically >g31 to Willaertia_magna|g31 - for some reason I dont really remeber, this is some strange exception.
grep ">g31" ./*_gappyout.fasta

#####################################################################################################
############################# CONCATENATING ######################################
#####################################################################################################
ls *_gappyout.fasta > alignments.txt

#replacing | with " "
for i in *_gappyout.fasta ; do cp "$i" ./"$i"_BACKUP.fasta ; done


# mkdir ../../../CONCAT_AND_GIANT_TREES_2_0/model_selection/
cp ./CONCAT* ../../../CONCAT_AND_GIANT_TREES_2_0/model_selection/


seqkit replace -p '(.+)' -r '{kv}' -k short_names.txt CONCAT.fa > renamedCONCAT.fa

./ElConcatenero.py -c -if fasta -of phylip -in renamedCONCAT.fa

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
#SBATCH -J cesgaPHYTESTconcat_pb_catgtr_ch1            # Job name
#SBATCH -o cesgaPHYTESTconcat_pb_catgtr_ch1_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e cesgaPHYTESTconcat_pb_catgtr_ch1_%j.e       # Name of stderr output file(%j expands to jobId)
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=1GB           # Cores per task requested
#SBATCH -t 100:00:00         # Run time (hh:mm:ss) - 2 days     

module load cesga/2020 gcc/system openmpi/4.1.1_ft3 phylobayes-mpi/1.9

cd /mnt/lustre/scratch/nlsas/home/csic/gfy/dri/dzavadska/pb_run_2_0/

srun --ntasks=64 pb_mpi -s -cat -gtr -d renamedCONCAT.phy cesgaPHYTESTconcat_pb_catgtr_ch3 
  

##############################################
##############################################
##############################################


sbatch -C ilk -t 100:00:00 --mem-per-cpu=1GB cesga_pb_w_renamedCONCAT_ch3.sh


##############################################
##############################################
##############################################


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

mpirun -np 64 pb_mpi cesgaPHYTESTconcat_pb_catgtr_ch5



##############################################
##############################################
##############################################



sbatch -C ilk -t 168:00:00 --mem-per-cpu=1GB cesga_pb_w_renamedCONCAT_RESUME_ch5.sh




###########################################
### tracecomp (if the chains converged)
#if needed, scp chains from another account 
cd /scratch/data2/dzavadska/COMPARATIVE_DATASET/CONCAT_AND_GIANT_TREES_2_0/pb_run_2_0/CESGA_run_2_0
scp csgfydri@ft3.cesga.es:/mnt/lustre/scratch/nlsas/home/csic/gfy/dri/dzavadska/pb_run_2_0_dri/* ./
scp csgfydaz@ft3.cesga.es:/home/csic/gfy/daz/dzavadska/pb_run_batch2/* ./


#on CESGA
module load cesga/2020 gcc/system openmpi/4.1.1_ft3 phylobayes-mpi/1.9

cd /mnt/lustre/scratch/nlsas/home/csic/gfy/dri/dzavadska/pb_run_2_0_dri/

bpcomp -x 4000 1 cesgaPHYTESTconcat_pb_catgtr_ch6 cesgaPHYTESTconcat_pb_catgtr_ch5 cesgaPHYTESTconcat_pb_catgtr_ch24 cesgaPHYTESTconcat_pb_catgtr_ch7 cesgaPHYTESTconcat_pb_catgtr_ch8 -o ./convergence/pb_2_0_5chains_Feb12




