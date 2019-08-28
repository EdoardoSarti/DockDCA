# 1) Gets the backmapping coverage from the corresponding folder in ddi_full_pdb/, and sorts the entries by coverage
# 2) Gets the resolution from the PDB (WARNING: it may not find the PDB or the resolution may not be there because of the experimental method used. In these cases, the entry is suppressed)
# 3) Writes a file with 

echo "#Pfam1_Pfam2 PDB_chain bkpcoverage1 bkpcoverage2 bkpcoverageTOT resolution" 
for x in `awk 'BEGIN{FS="_"} {print "ddi_"$2"_"$3}' pfam_pairs_version31_for_docking.txt`
do
	# echo ${x} 
	wc ddi_full_pdb/${x}/map1* | awk 'BEGIN{FS="_"} NF>7 {print $1, $7, substr($8,1,1)}' | awk -v x=${x} '{print $1, $5"_"$6}' > tmp1.txt
	wc ddi_full_pdb/${x}/map2* | awk 'BEGIN{FS="_"} NF>7 {print $1, $7, substr($8,1,1)}' | awk -v x=${x} '{print $1, $5"_"$6}' > tmp2.txt
	awk -v x=${x} 'BEGIN{while((getline < "tmp1.txt")>0) {a[$2]=$1}; while((getline < "tmp2.txt")>0) {b[$2]=$1}; for (i in a) {if (i in b) {print substr(x,5), i, a[i], b[i], a[i]+b[i]}}}' | sort -rnk5 > tmp3.txt
	N=`wc tmp3.txt | awk '{print $1}'`
	rm tmp4.txt
	rm tmp3bis.txt
	for ((i=1;i<=N;i++))
	do
		y=`awk -v i=${i} 'NR==i{print substr($2,1,4)}' tmp3.txt`
		grep "REMARK   2.*RESOLUTION" ../../../databases/PDB/${y}.pdb | awk '{print $4}' >> tmp4.txt
		EXST=`grep "REMARK   2.*RESOLUTION" ../../../databases/PDB/${y}.pdb | awk '{print $4}'`
		if [[ "${EXST}" != "" ]]
		then
			awk -v i=${i} 'NR==i' tmp3.txt >> tmp3bis.txt
		fi
	done
	paste tmp3bis.txt tmp4.txt -d " " > tmp5.txt
	# paste tmp3bis.txt tmp4.txt -d " " > check_${x}
	if [[ -s tmp5.txt ]]
	then
		awk 'NR==1{a=$5; m=$6; l=$0} NR>1&&$5==a&&$6<m{m=$6; l=$0} END {print l}' tmp5.txt | grep -v NOT
	fi
	STP=`grep 5nug tmp5.txt`
	if [[ "${STP}" != "" ]]
	then
		echo ${STP}
		exit
	fi
done
rm tmp*
