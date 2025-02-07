folder="/media/hard_disk/HDD_2_10TB/HaoJun/repeatlandscape/"
species=$1 # speices
term=$2
divs=(0 10 20 30 40)

original_bed=${species}.keeplong.gene.txt
TE_out=${species}.genome.fa.out

pre5k_outbed=${species}.Pre5k.bed
all_TE=${species}.TE_all.bed


### pre5k
if [ ! -e "${pre5k_outbed}" ]; then
	awk 'BEGIN{OFS="\t"} {
		if ($4 == "+") 
			if ($2 > 2) {
				{print $1, ($2 - 5000 > 0 ? $2 - 5000 : 1), $2-1} 
		} else if ($4 == "-")
			{print $1, $3+1, $3 + 5000}}' ${folder}${original_bed} | sort -k1,1 -k2,2n | bedtools merge > ${pre5k_outbed}
else
	printf "%4s   [Skip] generate %s\n" "" "${pre5k_outbed}"
fi

### All TE
for div in ${divs[@]}
do
	echo "${species}: div=${div}"
	if [[ "$term" == "All" || "$term" == "all" ]]; then
		awk -v term="$term" -v div="$div" 'BEGIN{OFS="\t"} {
			if (NR>3) {
				if ($2 >= div && $2 < (div + 10)) {
					if ($11 ~ /^DNA/ || $11 ~ /^LINE/ || $11 ~ /^LTR/){
						print $5, $6, $7, $11
					}
				}
			}
		}' ${folder}${TE_out} > ${species}.${term}.${div}.bed
	else
		awk -v term="$term" -v div="$div" 'BEGIN{OFS="\t"} {
			if (NR>3) {
				if ($2 >= div && $2 < (div + 10)) {
					if ($11 ~ "^" term){
						print $5, $6, $7, $11
					}
				}
			}
		}' ${folder}${TE_out} > ${species}.${term}.${div}.bed
	fi

	bedtools intersect -a ${pre5k_outbed} -b ${species}.${term}.${div}.bed -c > ${species}.count.${term}.${div}.bed
done