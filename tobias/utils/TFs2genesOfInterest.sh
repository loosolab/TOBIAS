#!/bin/bash

# author afust

# capability: 
# extracts all TF information for a target gene list

# parameters:
# target gene list: file with Ensembl gene ids in first column
# directory: directory containing files of interest (TFBS)
# pattern: file pattern ("*_overview.txt")
# output: output file


# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "	-g <string> []			file with Ensembl gene ids in first column"
	echo -e "	-d <string> []			directory containing files of interest (TFBS)"
	echo -e "	-p <string> []			file pattern ('*_overview.txt')"
	echo -e "   -o <string> []			output file"
	echo -e ""
	exit 1
}
export LC_ALL=C

genes=""
TFBS_dir=""
TFBS_pattern=""
output=""

while getopts ":g:d:p:o:" op; do
	case $op in
		g) genes=${OPTARG} ;;
		d) TFBS_dir=${OPTARG} ;;
		p) TFBS_pattern=${OPTARG} ;;
		o) output=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check parameters
if [[ -z $genes ]] || [[ -z $TFBS_dir ]] || [[ -z $TFBS_pattern ]] || [[ -z $output ]] ; then
	echo -e "input missing, please specify"
	usage
	exit -1
fi

# format gene list to expression used for grep (should not end with "|")
gene_list=$(awk 'BEGIN { ORS = "|" } { print }' ${genes})
if [ "${gene_list: -1}" == "|" ]; then
	gene_list=${gene_list::-1}
fi
# create header
f=$(find ${TFBS_dir} -name ${TFBS_pattern} | head -1)
head -1 ${f} > ${output}
# append all entries for all target genes to output
grep -h -E -r --include "${TFBS_pattern}" ${gene_list} ${TFBS_dir} >> ${output}