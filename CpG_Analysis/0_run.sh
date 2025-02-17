species=(RA	Pgle	BS	BP	PM	PMO	MC	RS	NM	Psem	SO	Stub)

term=$1

for sp in ${species[@]}
do
	printf "%4s => create5k_bed\n" "$sp"
	./1_create5k_bed.sh $sp $term

	printf "%4s => countGC\n" ""
	./2_creatPureInterval_and_CpG_cal.py $sp $term
done

echo ""

