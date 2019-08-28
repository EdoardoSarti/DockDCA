rm buried_sasa.txt
N=`wc best_coverage_and_resolution.txt | awk '{print $1}'`
for ((i=2;i<=N;i++))
do 
	PFs=`awk -v i=${i} 'NR==i{print $1}' best_coverage_and_resolution.txt`
	PDBC=`awk -v i=${i} 'NR==i{print $2}' best_coverage_and_resolution.txt`
	INT1=`ls ddi_full_pdb/*${PFs}*/map1_*_${PDBC}* | awk 'BEGIN{FS=":"}{print $NF; exit}'`
	INT2=`ls ddi_full_pdb/*${PFs}*/map2_*_${PDBC}* | awk 'BEGIN{FS=":"}{print $NF; exit}'`
	echo "${INT1} ${INT2}"
	ls ddi_full_pdb/*${PFs}*/map?_*_${PDBC}*
	PDB=`awk -v pdbc=${PDBC} 'BEGIN{print substr(pdbc, 1, 4)}'`
	PF1=`awk -v pfs=${PFs} 'BEGIN{print substr(pfs,1,7)}'`
	PF2=`awk -v pfs=${PFs} 'BEGIN{print substr(pfs,9,7)}'`
	awk -v int1=${INT1} -v pdbc=${PDBC} 'BEGIN {dashidx = index(int1, "-")
	                                            beg = substr(int1,1,dashidx-1)+0
	                                            term = substr(int1,dashidx+1)+0
	                                           }
	                                     $1=="ATOM"&&substr($0,22,1)==substr(pdbc,6,1)&&substr($0,23,4)+0>=beg&&substr($0,23,4)+0<=term
	                                    ' ../../../databases/PDB/${PDB}.pdb > interfaces/${PDBC}_${PF1}.pdb
	awk -v int2=${INT2} -v pdbc=${PDBC} 'BEGIN {dashidx = index(int2, "-")
	                                            beg = substr(int2,1,dashidx-1)+0
	                                            term = substr(int2,dashidx+1)+0
	                                           } 
	                                     $1=="ATOM"&&substr($0,22,1)==substr(pdbc,6,1)&&substr($0,23,4)+0>=beg&&substr($0,23,4)+0<=term
	                                    ' ../../../databases/PDB/${PDB}.pdb > interfaces/${PDBC}_${PF2}.pdb
	awk -v f1="interfaces/${PDBC}_${PF1}.pdb" -v f2="interfaces/${PDBC}_${PF2}.pdb" 'BEGIN {while((getline < f1)>0) {print substr($0,1,21)"A"substr($0,23)}; while((getline < f2)>0) {print substr($0,1,21)"B"substr($0,23)}}' > interfaces/${PDBC}_${PF1}_${PF2}_AB.pdb
	perl /home/sarti/Work/programs/ialign/bin/ialign.pl -w output interfaces/${PDBC}_${PF1}_${PF2}_AB.pdb AB interfaces/${PDBC}_${PF1}_${PF2}_AB.pdb AB -vmd cartoon -a 2 -dc 5 -do_self
	mv output/${PDBC}_${PF1}_${PF2}_ABAB_int.pdb interfaces/${PDBC}_${PF1}_${PF2}_AB_onlyINT.pdb
	vmd -dispdev text -e calculate_buried_sasa.tcl interfaces/${PDBC}_${PF1}_${PF2}_AB_onlyINT.pdb
	cat tmp_out_sasa.txt | awk -v foretxt="${PDBC}_${PF1}_${PF2}" '{printf "%s\t\t%10.2f\n", foretxt, $1}' >> buried_sasa.txt
done
sort -rnk2 buried_sasa.txt > interface_sasa.txt
rm buried_sasa.txt
rm -rf output/
rm *.vmd
rm tmp_out_sasa.txt

awk '$2>=300 {print $1}' interface_sasa.txt > list_cases_with_good_interface.txt
