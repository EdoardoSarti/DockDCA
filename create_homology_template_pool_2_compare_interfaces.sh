rm ialign_commands.txt
rm tm_scores.txt
rm conserved_interface.txt

echo "# TM-score %_orig_interface Min_P-value"

# Loop on the cases of good-interface Pfam domains
for x in `ls homology_template_pool/content_lists/contains_PF?????_wihtout_PF?????`
do
	W=`echo ${x} | awk 'BEGIN{FS="_"}{print $(NF-2)}'`
	WOUT=`echo ${x} | awk 'BEGIN{FS="_"}{print $(NF)}'`
	truefile=`ls interfaces/????_?_${W}_${WOUT}_AB_onlyINT.pdb`
	inv="no"
	if [[ "${truefile}" == "" ]]
	then
		truefile=`ls interfaces/????_?_${WOUT}_${W}_AB_onlyINT.pdb`
		inv="yes"
	fi
	truepdb=`echo ${truefile} | awk 'BEGIN{FS="/"} {print substr($2,1,4)}'`
	truech=`echo ${truefile} | awk 'BEGIN{FS="/"} {print substr($2,6,1)}'`

	# Loop on the candidate PDB structures for homology, containing the 1st PF (called WITH), and not the 2nd (called WITHOUT)
	for y in `cat ${x}`
	do
#		echo "${y} ${W} ${WOUT}"

		# dca2pdb finds the occurrences of WITH in the PDB y
		dca2pdb -pdb ${y} -b -v 31 > dca2pdb_out.txt
		awk -v w=${W} 'BEGIN {
		                      flag = 0
		                     }
		               $1=="Find" {
		                flag=0
		               }
		               flag==1&&index($1,w)!=0&&substr($1,1,2)=="PF" {
		                print
		               }
		               $1=="Backmap:" {
		                flag = 1
		               }' dca2pdb_out.txt > locations_with.txt

#		cat locations_with.txt
#		cat dca2pdb_out.txt
		N=`wc locations_with.txt | awk '{print $1}'`

		# If it has 0 locations, it means dca2pdb did not work properly
#		echo "${N} locs"
		if [[ "${N}" == "0" ]]
		then
			continue
		fi

		# Loop on the occurrences of the WITH domain in the PDB y
		for ((i=1;i<=N;i++))
		do
#			echo "loc ${i}"
			ch=`awk -v i=${i} 'NR==i{print $2}' locations_with.txt`
			intv=`awk -v i=${i} 'NR==i{print $3}' locations_with.txt`
			pfamacc=`awk -v i=${i} 'NR==i{print $1}' locations_with.txt`

			# 1) TM-SCORE
			# Isolate and write the current instance of the WITH domain in PDB y, in order to align it with the bound WITH domain
			awk -v ch=${ch} -v intv=${intv} 'BEGIN {dashidx = index(intv, "-")
			                                        beg = substr(intv,1,dashidx-1)+0
			                                        term = substr(intv,dashidx+1)+0
			                                       }
			                                 $1=="ATOM"&&substr($0,22,1)==ch&&substr($0,23,4)+0>=beg&&substr($0,23,4)+0<=term' ../../../databases/PDB/${y}.pdb |\
			awk '{printf "%s%5d%s\n", substr($0, 1, 6), NR, substr($0,12)}' > reference_pdbs/frtm_ref_${y}_${W}.pdb
		
			frf1="reference_pdbs/frtm_ref_${y}_${W}.pdb"
			frf2="interfaces/${truepdb}_${truech}_${W}.pdb"	

			# Align the current instance of the WITH domain in PDB y with the bound WITH domain
			# NOTE: Here, if the reference is too small, there can be warning flags
			# NOTE: The inversion in the arguments of the template and the target is done in order to have the correct TM-score (I THINK)
			/media/sarti/data/Work/programs/frtmalign/frtmalign ${frf2} ${frf1} -o superposition.pdb 2>&1 | grep "TM-score=" |\
			awk -v frf1=${frf1} -v frf2=${frf2} '{
			                                      for (i=1; i<=NF;i++) {
			                                          if (index($i, "TM-score")!=0) {
			                                              beg = index($i, "=")
			                                              term = index($i,",")
			                                              print frf1, frf2, substr($i, beg+1, term-beg-1)
			                                          }
			                                      }
			                                     }' > local_tm_score.txt
			cat local_tm_score.txt >> homology_template_pool/results/tm_scores.txt

			# 2) PRESENCE OF THE ORIGINAL INTERFACE	RESIDUES
			# Check how much of the interface between the bound WITH domain with the bound WITHOUT domain is reflected in the chosen instance of the WITH domain in PDB y
			awk -v w=${W} -v wout=${WOUT} -v pdb=${truepdb} -v ch=${truech} -v frf1=${frf1} -v frf2=${frf2} -v inv=${inv} 'BEGIN {
			                                        if (inv=="no") {
			                                            fx="interfaces/" pdb "_" ch "_" w "_" wout "_AB_onlyINT.pdb"
			                                        }
			                                        else
			                                        {
			                                            fx="interfaces/" pdb "_" ch "_" wout "_" w "_AB_onlyINT.pdb"
			                                        }
			                                        while((getline < fx)>0) {
			                                            if ($1=="ATOM"&&substr($0,22,1)=="A") {
			                                                interfaces[substr($0,23,4)+0]=-2
			                                            }
			                                        }
			                                 }
			                                 {
			                                  interfaces[substr($0,23,4)+0]++
			                                 } 
			                                 END {
			                                      for (i in interfaces) {
			                                          if (interfaces[i]==-1) {
			                                              n++
			                                              tot++
			                                          } 
			                                          else if (interfaces[i]==-2) {
			                                              tot++
			                                          }
			                                      }
			                                      printf "%s ORIGINAL_INTERFACE_RESIDUES_%s: %d %d %6.4f \n", frf1, wout, n, tot, n/tot
			                                 }' superposition.pdb > local_conserved_interface.txt
			rm superposition.pdb
			cat local_conserved_interface.txt >> homology_template_pool/results/conserved_interface.txt

			# P-VALUE OF THE INTERFACES
			# Looks for other instances of WITH in the same PDB chain
			rm others.txt
			awk -v ch=${ch} -v w=${W} -v pfamacc=${pfamacc} -v intv=${intv} 'BEGIN {
			                                                                        flag = 0
			                                                                        dashidx = index(intv, "-")
			                                                                        beg = substr(intv,1,dashidx-1)+0
			                                                                        term = substr(intv,dashidx+1)+0
			                                                                       }
			                                                                 $1=="Find" {
			                                                                  flag = 0
			                                                                 }
			                                                                 flag==1&&substr($1,1,2)=="PF"&&$2==ch&&index($1,pfamacc)==0 {
			                                                                  dashidxN = index($3, "-")
			                                                                  begN = substr($3,1,dashidxN-1)+0
			                                                                  termN = substr($3,dashidxN+1)+0
			                                                                  if (termN < beg || begN > term) {
			                                                                      print
			                                                                  }
			                                                                 } 
			                                                                 $1=="Backmap:" {
			                                                                  flag = 1
			                                                                 }' dca2pdb_out.txt > others.txt
#			cat others.txt

			if [[ -s others.txt ]]
			then
				# If there are other occurrences of WITH or other PF domains in the same PDB chain,
				# we have to separate the wanted occurrence from the rest of the chain, assigning it the chain name "0". 
				# The chain now will contain a hole: the one created by singling out WITH.
				awk -v ch=${ch} -v intv=${intv} 'BEGIN {dashidx = index(intv, "-")
				                                        beg = substr(intv,1,dashidx-1)+0
				                                        term = substr(intv,dashidx+1)+0
				                                       }
				                                 {
				                                  if ($1=="ATOM"&&substr($0,22,1)==ch&&substr($0,23,4)+0>=beg&&substr($0,23,4)+0<=term) {
				                                      n++
				                                      l[n] = substr($0,1,21)"0"substr($0,23)
				                                  }
				                                  else if ($1=="ATOM")
				                                  {
				                                      print
				                                  }
				                                 }
				                                 END {
				                                  for (i=1;i<=n;i++) {
				                                      print l[i]
				                                  }
				                                 }' ../../../databases/PDB/${y}.pdb |\
				awk '{printf "%s%5d%s\n", substr($0, 1, 6), NR, substr($0,12)}' > reference_pdbs/ref_${y}.pdb

				# And then we take note of the interfaces to check
				awk 'BEGIN {
				            flag = 0
				           }
				$1=="Find" {
				 flag = 0
				}
				flag==1&&substr($1,1,2)=="PF" {
				 a[$2] = 1
				}
				$1=="Backmap:" {
				 flag = 1
				}
				END {
				     for (i in a) {
				         print "0"i
				     }
				    }' dca2pdb_out.txt > interfaces_to_check.txt

			else
				# Otherwise, we just take note of the interfaces to check
				awk -v pfamacc=${pfamacc} -v ch=${ch} 'BEGIN {
				                                              flag = 0
				                                             }
				                                       $1=="Find" {
				                                        flag = 0
				                                       }
				                                       flag==1&&substr($1,1,2)=="PF"&&index($1,pfamacc)==0 {
				                                        a[$2] = 1
				                                       }
				                                       $1=="Backmap:" {
				                                        flag = 1
				                                       } 
				                                       END {
				                                        for (i in a) {
				                                            print ch""i
				                                        }
				                                       }' dca2pdb_out.txt > interfaces_to_check.txt

				cp ../../../databases/PDB/${y}.pdb reference_pdbs/ref_${y}.pdb
			fi

			rm others.txt

			for interf in `cat interfaces_to_check.txt`
			do
				echo "perl /home/sarti/Work/programs/ialign/bin/ialign.pl -w output ${truefile} AB reference_pdbs/ref_${y}.pdb ${interf} -vmd cartoon -a 2 -dc 5" > local_ialign_commands.txt
			done
			rm interfaces_to_check.txt
			cat local_ialign_commands.txt >> ialign_commands.txt

			Ncomm=`wc local_ialign_commands.txt | awk '{print $1}'`
			for ((j=1;j<=Ncomm;j++))
			do 
				awk -v j=${j} 'NR==j' local_ialign_commands.txt
				`awk -v j=${j} 'NR==j' local_ialign_commands.txt` | grep "P-value"
			done > pre_p_values.txt
			rm *.vmd

			awk '$1=="perl"{x=$7" "$8" "$5" "$6} $1=="IS-score"{printf "%s %10.6f\n", x, substr($6,1,length($6)-1)}' pre_p_values.txt > local_p_values.txt
			rm pre_p_values.txt
			cat local_p_values.txt >> homology_template_pool/results/p_values.txt

			TMSCORE=`awk '{print $NF}' local_tm_score.txt`
			ORIGINT=`awk '{print $NF}' local_conserved_interface.txt`
			MINPVALUE=`awk '{print $NF}' local_p_values.txt | sort | head -n1`

			if [[ "${MINPVALUE}" == "" ]]
			then
				MINPVALUE="ERR"
			fi

			echo "${W} without ${WOUT} in ${y} ${ch} (${pfamacc}): ${TMSCORE} ${ORIGINT} ${MINPVALUE}"
			rm local_*
		done
		rm dca2pdb_out.txt locations_with.txt
	done
done
