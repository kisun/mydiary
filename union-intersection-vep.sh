print the most common obs in a tabular file
awk '{print $1}' deleterious_sorted.txt | sort | uniq -c | sort -rn | awk '{if ($1 >=10) print ($1 "\t" $2)}'| wc -l


Extract all deleterious SNPs from the vep prediction
grep "deleterious" vep_raw_filter_715_5_NonSyn.txt|sort > deleterious.txt (we could already sort the file using sort)

Extract all deleterious SNPs associated with Finnsheep
egrep -r 'IND=21J|IND=23C|IND=35C|IND=39A|IND=39E|IND=39J|IND=5B|IND=801|IND=41B|IND=42B|IND=45A|IND=46F|IND=974' deleterious_sort.txt > deleterious_FIN.txt

Extract all deleterious SNPS associated with Texel
egrep -r 'IND=13B|IND=22A|IND=379|IND=43A|IND=4H|IND=4P|IND=787|IND=9K|IND=11E|IND=26B|IND=27A|IND=3B|IND=44A' deleterious_sort.txt > deleterious_TEX.txt

Extract all deleterious SNPs associated with F1
egrep 'IND=24A|IND=31A|IND=33A|IND=7G|IND=7M|IND=19A|IND=29A|IND=37C5|IND=37C6|IND=37C7|IND=4208|IND=4563|IND=4590' deleterious_sort.txt | wc -l

Extract all deleterious SNPs (only IDs from Finnsheep, texel and F1) if they are present in at least 2 samples and 5 samples
awk '{print $1}' deleterious_FIN.txt | sort | uniq -c | sort -rn |awk '{if ($1 >=2) print ($2)}'> deleterious_FIN_geq2.txt
awk '{print $1}' deleterious_FIN.txt | sort | uniq -c | sort -rn |awk '{if ($1 >=5) print ($2)}'> deleterious_FIN_geq5.txt

and so on

comm -12 deleterious_FIN_geq5.txt deleterious_TEX_geq5.txt > deleterious_FIN_TEX_common.txt #Common to FIN and TEX (including SNPS that are common to all i.e intersection of three breeds)
comm -23 deleterious_FIN_geq5.txt deleterious_TEX_geq5.txt > deleterious_FINnotTEX.txt #SNPS in FIN that are not in TEX
comm -23 deleterious_TEX_geq5.txt deleterious_FIN_geq5.txt > deleterious_TEXnotFIN.txt #SNPS in TEX that are not in FIN

and so on..

comm -23 deleterious_FXTnotFIN.txt deleterious_TEX_geq5.txt > deleterious_onlyFXT.txt #deleterious SNPS that are unique to FXT
comm -23 deleterious_TEXnotFXT.txt deleterious_FIN_geq5.txt > deleterious_onlyTEX.txt #deleterious only in TEX
comm -23 deleterious_FINnotFXT.txt deleterious_TEX_geq5.txt > deleterious_onlyFIN.txt #deleterious only in FIN

#Identify common SNPS in two breeds excluding the SNPS that are common to all. i.e SNPS that are common to only given two breeds (i.e excluding intersection of three breeds)
comm -12 deleterious_FXTnotFIN.txt deleterious_TEXnotFIN.txt > deleterious_comm_TEX_FXT.txt #common to TEX and FXT (excluding SNPS that are common to all)
comm -12 deleterious_FINnotTEX.txt deleterious_FXTnotTEX.txt > deleterious_comm_FIN_FXT.txt#
comm -12 deleterious_FINnotFXT.txt deleterious_TEXnotFXT.txt > deleterious_comm_FIN_TEX.txt

Identify common SNPS in all breeds
comm -23 deleterious_FIN_TEX_common.txt deleterious_comm_FIN_TEX.txt > deleterious_comm_FIN_TEX_FXT.txt

Now extract detail SNP information using unique SNP ids from above analysis.
awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' deleterious_FINnotTEX.txt deleterious_FIN.txt > deleterious_FINnotTEX_detail.txt #extract all SNPs detail for the SNPS that are in FIN but not in TEX.
and so on..

Second Round
Extract unique IDS from FIN, TEX and FXT deleterious SNPS
awk '{print $1}' deleterious_FIN.txt| sort | uniq -c|awk '{print $2}'> deleterious_FIN_id.txt 
awk '{print $1}' deleterious_TEX.txt| sort | uniq -c|awk '{print $2}'> deleterious_TEX_id.txt 
awk '{print $1}' deleterious_FXT.txt| sort | uniq -c|awk '{print $2}'> deleterious_FXT_id.txt 

FIND SNPS in one breed that are not in another breed
comm -23 deleterious_FIN_id.txt deleterious_TEX_id.txt > deleterious_FINnotTEX.txt
comm -23 deleterious_TEX_id.txt deleterious_FIN_id.txt > deleterious_TEXnotFIN.txt
comm -23 deleterious_FIN_id.txt deleterious_FXT_id.txt > deleterious_FINnotFXT.txt
comm -23 deleterious_FXT_id.txt deleterious_FIN_id.txt > deleterious_FXTnotFIN.txt
comm -23 deleterious_TEX_id.txt deleterious_FXT_id.txt > deleterious_TEXnotFXT.txt
comm -23 deleterious_FXT_id.txt deleterious_TEX_id.txt > deleterious_FXTnotTEX.txt

FIND common SNPS in two breeds
comm -12 deleterious_FIN_id.txt deleterious_TEX_id.txt > deleterious_comm_FIN_TEX.txt
comm -12 deleterious_FIN_id.txt deleterious_FXT_id.txt > deleterious_comm_FIN_FXT.txt
comm -12 deleterious_TEX_id.txt deleterious_FXT_id.txt > deleterious_comm_TEX_FXT.txt

FIND common SNPS in three breeds
comm -12 deleterious_comm_FIN_TEX.txt deleterious_FXT_id.txt > deleterious_comm_FIN_TEX_FXT.txt 

SNPs found only in given breed
comm -23 deleterious_FINnotTEX.txt deleterious_FXT_id.txt > deleterious_onlyFIN.txt
comm -23 deleterious_TEXnotFIN.txt deleterious_FXT_id.txt > deleterious_onlyTEX.txt
comm -23 deleterious_FXTnotFIN.txt deleterious_TEX_id.txt > deleterious_onlyFXT.txt
