# Using-viromes-to-predict-novel-immune-proteins-in-non-model-organisms
The following commands were used to process and analyze viromes from manuscript:

Quistad SD, Lim YW, Gueiros Z. Silva G., Nelson CE, Haas AF, Wegley Kelly L, Edwards RA, and Rohwer FL (2016). Using viromes to predict novel immune proteins in non-model organisms. Proceedings of the Royal Society B.

#####################################
VIROME PREPROCESSING
#####################################
Deconseq was run against all viromes with the following parameters
#####################################

$ perl deconseq.pl -f Acr1Fe_P.fasta -dbs AdigGenMito -i 90 -c 90 -S 10000000 -z 1 -T 30 -id Acr1Fe_PP

#####################################
Prinseq was run against all viromes to remove exact duplicates, change fasta headers to include biological sample, change fasta file names to _F (F=Final for downstream analysis)
#####################################

$ perl prinseq-lite.pl -fasta Acr1LT_PP_clean.fa -derep 1

$ perl prinseq-lite.pl -fasta Acr1LT_F.fa -seq_ID Acr1LT_

########################################
All viruses were combined into single fasta using the cat command (CoralVir.fa) 1,048,627 sequences total representing 14 individual corals across 4 species.
########################################



#####################################
The following commands were used to perform VIROME AND CDD ANALYSIS
########################################
rpsblast coral proteins vs CDD database
########################################

$ nohup rpsblast -i adi_prot.fa -d /home/squistad/CDD/CDDdb/Cdd -o adi_protVsCDD.m8 -e 0.01 -m 8 &

########################################
tBLASTn of all coral viruses generated here (1048627 sequences) 
########################################

$ nohup tblastn -db CoralVir.fa -query adi_prot.fa -evalue 1E-04 -seg yes -comp_based_stats D -outfmt 6 -out CoralVirVsAdiProt.m8 &

########################################
Extract out top hit of each viral sequence to one coral protein
########################################

$ sort -k2,2 -rk12 CoralVirVsAdiProt.m8 | awk 'n!=$2{print;n=$2}' > CoralVirVsAdiProt.m8_uniq

########################################
Count # viral hits to coral proteins
########################################

$ cut -f1 CoralVirVsAdiProt.m8_uniq | sort | uniq -c | sort -nr > CoralVirVsAdiProt.m8_uniq_counts

########################################
Create two columns in count file
########################################

$ tr ' ' '\t' < CoralVirVsAdiProt.m8_uniq_counts > CoralVirVsAdiProt.m8_uniq_counts_tab

########################################
Extract out proteins with viral hits
########################################

$ cut -f1 CoralVirVsAdiProt.m8_uniq | sort | awk 'n!=$1{print;n=$1}' > CoralVirVsAdiProt.m8_uniq_f1

$ perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' CoralVirVsAdiProt.m8_uniq_f1 adi_prot.fa > CoralVirVsAdiProt.fa

########################################
Blast proteins with viral hits against the CDD
########################################

$ nohup rpsblast -i CoralVirVsAdiProt.fa -d /home/squistad/CDD/CDDdb/Cdd -o CoralVirVsAdiProt_CDD -e 0.01 -m 8 &

########################################
Blast proteins against themselves
########################################

$ nohup blastp -db CoralVirVsAdiProt.fa -query CoralVirVsAdiProt.fa -evalue 1E-10 -seg yes -comp_based_stats D -outfmt 6 -out CoralVirProtvsCoralVirProt.m8 &

########################################
Extract out top hits
########################################

$ sort -k2,2 -rk12 CoralVirProtvsCoralVirProt.m8 | awk 'n!=$2{print;n=$2}' > CoralVirProtvsCoralVirProt.m8_uniq



########################################
VIRUS GENOMIC ANALYSIS (HUMAN)
########################################
Split human virus genomes into 200bp chunks in single fasta file (Human herpesvirus 1, Human herpesvirus 2 strain HG52, Human_herpesvirus_3_uid15198, Human_herpesvirus_4_uid14413, Human_herpesvirus_4_uid14413, Human_herpesvirus_6A_uid14462,  Human_herpesvirus_6B_uid14422, Human_herpesvirus_7_uid14625, Human_herpesvirus_7_uid14625, Human_adenovirus_A_uid14517, Human_adenovirus_A_uid14517, Human_adenovirus_C_uid14518, Human_adenovirus_D_uid14535, Human_adenovirus_E_uid15152, Human_circovirus_VS6600022_uid257742, Human_papillomavirus___1_uid15491, Human_papillomavirus___2_uid15512
########################################

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirAGen.fa > AdenoVirAGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirBGen.fa > AdenoVirBGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirCGen.fa > AdenoVirCGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirDGen.fa > AdenoVirDGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirEGen.fa > AdenoVirEGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' AdenoVirFGen.fa > AdenoVirFGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir1Gen.fa > HerpVir1Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir2Gen.fa > HerpVir2Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir3Gen.fa > HerpVir3Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir4Gen.fa > HerpVir4Gen.fa_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir5Gen.fa > HerpVir5Gen.fa_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir6AGen.fa > HerpVir6AGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir6BGen.fa > HerpVir6BGen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir7Gen.fa > HerpVir7Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HerpVir8Gen.fa > HerpVir8Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HPV1Gen.fa > HPV1Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' HPV2Gen.fa > HPV2Gen_split.fa

$ perl -ne 'BEGIN{$/=">"}{s/(.*)//;$n=$1;s/\n//g;$i=0;s/(.{1,200})/printf(">%s_%05d\n%s\n",$n,++$i,$1);/ge;}' CircVirGen.fa > CircVirGen_split.fa

########################################
Rename fasta headers
########################################

$ perl prinseq-lite.pl -fasta AdenoVirAGen_split.fa -seq_ID AdenoVirA -out_good AdenoVirAGen_split_newids

$ perl prinseq-lite.pl -fasta AdenoVirBGen_split.fa -seq_ID AdenoVirB -out_good AdenoVirBGen_split_newids

$ perl prinseq-lite.pl -fasta AdenoVirCGen_split.fa -seq_ID AdenoVirC -out_good AAdenoVirCGen_split_newids

$ perl prinseq-lite.pl -fasta AdenoVirDGen_split.fa -seq_ID AdenoVirD -out_good AdenoVirDGen_split_newids

$ perl prinseq-lite.pl -fasta AdenoVirEGen_split.fa -seq_ID AdenoVirE -out_good AdenoVirEGen_split_newids

$ perl prinseq-lite.pl -fasta AdenoVirFGen_split.fa -seq_ID AdenoVirF -out_good AdenoVirFGen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir1Gen_split.fa -seq_ID HerpVir1 -out_good HerpVir1Gen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir2Gen_split.fa -seq_ID HerpVir2 -out_good HerpVir2Gen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir3Gen_split.fa -seq_ID HerpVir3 -out_good HerpVir3GenGen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir4Gen.fa_split.fa -seq_ID HerpVir4 -out_good HerpVir4Gen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir5Gen.fa_split.fa -seq_ID HerpVir4 -out_good HerpVir5Gen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir6AGen_split.fa -seq_ID HerpVir6A -out_good HerpVir6AGen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir6BGen_split.fa -seq_ID HerpVir6B -out_good HerpVir6BGen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir7Gen_split.fa -seq_ID HerpVir7 -out_good HerpVir7Gen_split_newids

$ perl prinseq-lite.pl -fasta HerpVir8Gen_split.fa -seq_ID HerpVir8 -out_good HerpVir8Gen_split_newids

$ perl prinseq-lite.pl -fasta HPV1Gen_split.fa -seq_ID HPV1 -out_good HPV1Gen_split_newids

$ perl prinseq-lite.pl -fasta HPV2Gen_split.fa -seq_ID HPV2 -out_good HPV2Gen_split_newids

$ perl prinseq-lite.pl -fasta CircVirGen_split.fa -seq_ID CircVir -out_good CircVirGen_split_newids

########################################
Combine into single fasta
########################################

$ cat AdenoVirAGen_split_newids.fasta AdenoVirBGen_split_newids.fasta AAdenoVirCGen_split_newids.fasta AdenoVirDGen_split_newids.fasta AdenoVirEGen_split_newids.fasta AdenoVirFGen_split_newids.fasta CircVirGen_split_newids.fasta HerpVir1Gen_split_newids.fasta HerpVir2Gen_split_newids.fasta HerpVir3GenGen_split_newids.fasta HerpVir4Gen_split_newids.fasta HerpVir5Gen_split_newids.fasta  HerpVir6AGen_split_newids.fasta HerpVir6BGen_split_newids.fasta HerpVir7Gen_split_newids.fasta HerpVir8Gen_split_newids.fasta HPV1Gen_split_newids.fasta HPV2Gen_split_newids.fasta > HumVirAll_200bp.fa

########################################
tblastn HumVirAll_200bp.fa vs human proteome
########################################

$ nohup tblastn -db HumVirAll_200bp.fa -query H.sapien_prot.fa -evalue 1E-04 -seg yes -comp_based_stats D -outfmt 6 -out HumVirAll_200bpVsProt.m6 &

########################################
Extract out top hits
########################################

$ sort -k2,2 -rk12 HumVirAll_200bpVsProt.m6 | awk 'n!=$2{print;n=$2}' > HumVirAll_200bpVsProt.m6_uniq

########################################
Count # viral hits to human proteins
########################################

$ cut -f1 HumVirAll_200bpVsProt.m6_uniq | sort | uniq -c | sort -nr > HumVirAll_200bpVsProt.m6_uniq_counts

########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ cut -f1 HumVirAll_200bpVsProt.m6_uniq | sort | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > HumVirAll_200bpVsProt_sort_nonumbers

########################################
Remove lines from one file which are also in another file (removing hits to exclusion_protein_database)
########################################

$ awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 HumMetabProt_UniprotIDs f=2 HumVirAll_200bpVsProt_UniprotIDs > HumVirAll_200bpVsProt_UniprotIDs_filtered

########################################
Count occurrences of filtered hits
########################################

$ sort HumVirAll_200bpVsProt_AllHits_Filtered | uniq -c | sort -nr > HumVirAll_200bpVsProt_Filtered_counts

########################################
HUMAN VIRAL METAGENOMIC ANALYSIS (HUMAN)
########################################
Final proteins for removal from BLAST results
########################################

$ cat HumVirVsProt_BLAST_prot_rem_UniprotIDs Cytoskeletal_structural_UnitProtIDs  Mito_UniProtIDs HumMetabProt_UniprotIDs Ribonuclear_UniProtIDs > UniProtIDs_Rem

$ sort UniProtIDs_Rem | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > UniProtIDs_Rem_uniq

$ UniProtIDs_Rem_uniq = 8946 proteins total

$ mv UniProtIDs_Rem_uniq Human_UniProtIDs_Rem_uniq

########################################
Remove lines from one file which are also in another file (removing hits to exclusion_protein_database).
########################################

$ awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 UniProtIDs_Rem_uniq f=2 HumVirVsProt.m6_uniq_f1_tab_f1 > HumVirVsProt_UniprotIDs_Filtered

########################################
Count occurrences of filtered hits
########################################

$ sort HumVirVsProt_UniprotIDs_Filtered | uniq -c | sort -nr > HumVirvsProt_Filtered_counts

########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ sort HumVirVsProt_UniprotIDs_Filtered  | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > HumVirVsProt_UniprotIDs_Filtered_counts_nonumbers

########################################
Count occurrences of filtered hits
########################################

$ sort HumVirVsProt_UniprotIDs_Filtered | uniq -c | sort -nr > HumVirVsProt_UniprotIDs_Filtered_counts

########################################
BLASTp human proteins for removal vs coral proteome
########################################

$ nohup blastp -db Human_prot_removal.fa -query adi_prot.fa -evalue 1E-05 -seg yes -comp_based_stats D -outfmt 6 -out Human_prot_remVSCoralProt.m6 

########################################
Extract out top hits 
########################################

$ sort -k2,2 -rk12 Human_prot_remVSCoralProt.m6 | awk 'n!=$2{print;n=$2}' > Human_prot_remVSCoralProt_uniq


########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ sort CoralVirVsAdiProt.m8_uniq_f1  | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > CoralVirVsAdiProt.m8_uniq_f1_sort_nonumbers

########################################
Remove hits to database from coral hits
########################################

$ awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 Human_prot_remVSCoralProt_uniq_coralIDs f=2 CoralVirVsAdiProt.m8_uniq_f1 > CoralVirVsAdiProt_IDs_Filtered

########################################
Extract out top hits
########################################

$ sort -k2,2 -rk12 CoralVirVsAdiProt_IDs_Filtered | awk 'n!=$2{print;n=$2}' > CoralVirVsAdiProt_IDs_Filtered

########################################
Count occurrences of filtered hits
########################################

$ sort CoralVirVsAdiProt_IDs_Filtered  | uniq -c | sort -nr > CoralVirVsAdiProt_IDs_Filtered_counts

########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ sort CoralVirVsAdiProt_IDs_Filtered | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > CoralVirVsAdiProt_IDs_Filtered_counts_nonumbers

########################################
PROCESSING CORAL VIROME HITS UNFILTERED
########################################

$ sort CoralVirVSAdProt.m6_uniq_f1  | uniq -c | sort -nr > CoralVirVSAdProt.m6_uniq_f1_counts

########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ sort CoralVirVSAdProt.m6_uniq_f1 | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > CoralVirVSAdProt.m6_uniq_f1_counts_nonumbers

########################################
Blastp of coral proteins with at least 4 hits against the human proteome
########################################

$ nohup blastp -db /home/squistad/ViralPopImmPred/BLAST/H.sapien_prot.fa -query CoralVirVSAdProt.fa -evalue 1E-05 -seg yes -comp_based_stats D -outfmt 6 -out CoralVirVSAdProtVSHumanProt.m6 &

########################################
Extract out protein hits
########################################

$ cut -f1 CoralVirVSAdProtVSHumanProt.m6 > CoralVirVSAdProtVSHumanProt.m6_f1

$ sort CoralVirVSAdProtVSHumanProt.m6_f1 | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > CoralVirVSAdProtVSHumanProt.m6_f1_uniq

########################################
Remove proteins that have human hits
########################################

$ awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 CoralVirVSAdProtVSHumanProt.m6_f1_uniq f=2 CoralVirVSAdProt_IDS > CoralVirVSAdProtVSHumanProt_no_human_hits

########################################
Extract out sequence files from a fasta based on separate file that contains a list of all the things you want to extract in a single column
########################################

$ perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' CoralVirVSAdProtVSHumanProt_no_human_hits CoralVirVSAdProt.fa > CoralVirVSAdProtVSHumanProt_no_human_hits.fa

########################################
Blast check
########################################

$ nohup blastp -db /home/squistad/ViralPopImmPred/BLAST/H.sapien_prot.fa -query CoralVirVSAdProtVSHumanProt_no_human_hits.fa -evalue 1E-05 -seg yes -comp_based_stats D -outfmt 6 -out CoralVirVSAdProtVSHumanProt_no_human_hitsVSHumanProt.m6 &

########################################
Extract out top hits
########################################

$ sort -k2,2 -rk12 CoralVirVSAdProtVSHumanProt.m6 | awk 'n!=$2{print;n=$2}' > CoralVirVSAdProtVSHumanProt.m6_uniq

########################################
Extract out rows based on sequence ID file (for example extract out rows from m6 blast file from another file that contains all the protein IDs you are interested in)
########################################

$ grep -f CoralVirVSAdProtVSHumanProt_no_human_hits CoralVirVsAdiProt.m8_uniq > CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full

$ cut -f1 CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full > CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full_f1

$ sort CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full_f1  | uniq -c | sort -nr > CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full_f1_counts

########################################
Sort file which contains 1 column and multiple occurrences of ID by most abundant to least abundant and print without actual number displayed
########################################

$ sort CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full_f1 | uniq -c | sort -rn | sed -E 's/^ *[0-9]+ //g' > CoralVirVSAdProtVSHumanProt_no_human_hits_m6_full_f1_counts_nonumbers

########################################
rpsblast coral proteins that DO NOT have human hits vs CDD database
########################################

$ nohup rpsblast -i CoralVirVSAdProtVSHumanProt_no_human_hits.fa -d /home/squistad/CDD/CDDdb/Cdd -o CDD_CoralVirVSAdProtVSHumanProt_no_human_hits.m8 -e 0.01 -m 8 &
