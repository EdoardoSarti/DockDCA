for x in `cat list_cases_with_good_interface.txt`
do
	echo ${x} 
	PF1=`echo ${x} | awk 'BEGIN{FS="_"}{print $3}'`
	PF2=`echo ${x} | awk 'BEGIN{FS="_"}{print $4}'`
	zgrep ${PF1} ../../pfam_interactions/db/external_resources/Pfam_31/database_files/pdbmap.gz | awk '{print tolower(substr($1,1,4))}' | sort | uniq > homology_template_pool/content_lists/contains_${PF1}
	zgrep ${PF2} ../../pfam_interactions/db/external_resources/Pfam_31/database_files/pdbmap.gz | awk '{print tolower(substr($1,1,4))}' | sort | uniq > homology_template_pool/content_lists/contains_${PF2}
	grep -v -f homology_template_pool/content_lists/contains_${PF2} homology_template_pool/content_lists/contains_${PF1} > homology_template_pool/content_lists/contains_${PF1}_wihtout_${PF2}
	grep -v -f homology_template_pool/content_lists/contains_${PF1} homology_template_pool/content_lists/contains_${PF2} > homology_template_pool/content_lists/contains_${PF2}_wihtout_${PF1}
	for y in `cat homology_template_pool/content_lists/contains_${PF1}_wihtout_${PF2}`
	do
		rm homology_template_pool/content_lists/only_${PF1}_in_pdb
		yy=`echo "${y}" | awk '{print toupper(substr($1,1,4))}'`
		N1=`zgrep ${yy} ../../pfam_interactions/db/external_resources/Pfam_31/database_files/pdbmap.gz | wc | awk '{print $1}'`
		N2=`zgrep ${PF1} ../../pfam_interactions/db/external_resources/Pfam_31/database_files/pdbmap.gz | grep ${yy} | wc | awk '{print $1}'`
		if [[ "${N1}" == "${N2}" ]]
		then
			echo "${y}" >> homology_template_pool/content_lists/only_${PF1}_in_pdb
		fi
	done
done
