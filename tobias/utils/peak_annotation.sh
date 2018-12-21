#!/bin/bash

# author afust

# capabilities:
# peak annotation

# requirements: 
# UROPA
# py2
# R

usage() {
	echo -e ""
	echo -e "Usage: $0 [options]"
	echo -e "Mandatory"
	echo -e "	-i|--bed <string> []			input bed file for annotation"
	echo -e "	-o|--annotation <string> []		peak annotation file"
	echo -e "	-g|--gtf <string> []			gtf file for annotation, check filter attributes"
	echo -e "Optional"
	echo -e "	-p|--header <string> [bed_header.txt]	peak annotation header file"
	echo -e "	-l|--log <string> [input.log]		log file"
	echo -e "	-c|--cores <int> [1]			number of cores"
	#echo -e "	-j|--universe <string> [] gene universe, can be used for enrichment"
	echo -e "For UROPA params check: https://uropa-manual.readthedocs.io/config.html"
	echo -e "	-t|--feature <string> [gene]		UROPA param"
	echo -e "	-r|--feature_anchor <string> [start]	UROPA param"
	echo -e "	-e|--distance <string> [[10000,1000]]	UROPA param"
	echo -e "	-a|--filter_attribute <string> []	UROPA param"
	echo -e "	-w|--attribute_value <string> []	UROPA param"
	echo -e "	-n|--show_attribute <string> [gene_name,gene_id,gene_biotype]	UROPA param"
	exit 1
}
export LC_ALL=C

# defaults
cores=1
feature="gene"
feature_anchor="start"
distance=[10000,1000]
#filter_attribute="gene_biotype"
#attribute_value="protein_coding"
show_attribute='"gene_biotype", "gene_name", "gene_id"'

# no option given
if [ $# -eq 0 ]; then
	usage
	exit -1
fi

# Process command line options
progname=${0##*/}
shortopts="i:o:j:p:l:c:g:t:r:e:a:w:n:h"
longopts="bed:,annotation:,universe:,header:,log:,cores:,gtf:,feature:,feature_anchor:,distance:,filter_attribute:,attribute_value:,show_attribute:"
opts=$(getopt -s bash --options "$shortopts" --longoptions "$longopts" --name "$progname" -- "$@")

eval set -- ${opts}
while true; do
    case "$1" in
        ( --                    ) break ;;
        ( -h|--help             ) usage ;;
        ( -i|--bed              ) bed=$2; shift ;;
        ( -o|--annotation       ) annotation=$2; shift ;;
        ( -j|--universe          ) universe=$2; shift ;;
        ( -p|--header        	) header=$2; shift ;;
        ( -l|--log              ) log=$2; shift ;;
        ( -c|--cores            ) cores=$2; shift ;;
        ( -g|--gtf              ) gtf=$2; shift ;;
        ( -t|--feature          ) feature=$2; shift ;;
        ( -r|--feature_anchor   ) feature_anchor=$2; shift ;;
	    ( -e|--distance         ) distance=$2; shift ;;
        ( -a|--filter_attribute	) filter_attribute=$2; shift ;;
        ( -w|--attribute_value  ) attribute_value=$2; shift ;;
        ( -n|--show_attribute   ) show_attribute=$2; shift ;;
    esac
    shift
done
# process params
if [[ -z $bed ]] || [[ -z $annotation ]] || [[ -z $gtf ]]; then
	echo -e "bed file or gtf for annotation or output missing: mandatory!" 
	usage
	exit -1
fi
if [[ -z $header ]]; then
	tmp=$(cut -f1 -d'.' $(basename ${annotation}))'_header.txt'
	header=$(dirname ${annotation})'/'$tmp
fi
if [[ -z $log ]]; then
	log=$(basename $bed .bed)'.log'
fi

if [[ ${show_attribute:0:1} != '"' ]]; then
	show_attribute=$(echo $show_attribute | sed -r 's/[,]/","/g' | sed -e 's/.*/"&"/')
fi

prefix="$(dirname $annotation)/$(basename $annotation .bed)"
json="${prefix}.json"


echo -e "Started peak annotation:" 2>&1 | tee -a $log
# create json config for uropa annotation
if [[ ! -z $filter_attribute ]] && [[ ! -z $attribute_value ]]; then
	cat > $json << EOL
{
	"queries":[
	{"feature":"$feature", "feature.anchor":"$feature_anchor", "distance":$distance, "filter.attribute":"$filter_attribute", "attribute.value":"$attribute_value", "show.attributes":[$show_attribute]},
	{"feature":"$feature", "feature.anchor":"$feature_anchor", "distance":$distance, "show.attributes":[$show_attribute]}
	],
	"priority": "True",
	"gtf": "$gtf",
	"bed": "$bed"
}
EOL
else
	cat > $json << EOL
{
	"queries":[
	{"feature":"$feature", "feature.anchor":"$feature_anchor", "distance":$distance, "show.attributes":[$show_attribute]}
	],
	"gtf": "$gtf",
	"bed": "$bed"
}
EOL
fi

# uropa annotation
echo -e "uropa -i $json -l $log -t $cores -p $prefix -d -s" 2>&1 | tee -a $log
uropa -i $json -l $log -t $cores -p $prefix -d -s 2>&1 | tee -a $log


# build peaks with annotation
# create final header of all merged annotated
attributes=$(($(head -1 ${prefix}_finalhits.txt | awk '{print NF}')-1))
echo -e "attributes: $attributes"
`head -1 ${prefix}_finalhits.txt | awk 'BEGIN { OFS = "\n" } ; { print $2,$3,$5,$1,$6,$7,$8,$9,$10,$11 }' - > $header`
`head -1 ${prefix}_finalhits.txt | cut -f15-${attributes} - | sed -e 's/\t/\n/g' >> $header`
# extract columns of interest
`awk 'BEGIN { OFS = "\t" } ; { print $2,$3,$5,$1,$6,$7,$8,$9,$10,$11 }' ${prefix}_finalhits.txt | sed -e 1d |  sed -e 's/\t\t/\t/g' > $annotation'.tmp'` 	
cut -f15-${attributes} ${prefix}_finalhits.txt | sed -e 1d | paste -d"\t" $annotation'.tmp' - > $annotation
sort -k1,1 -k2,2n $annotation > $annotation'.tmp'
mv $annotation'.tmp' $annotation

# add in case of go enrichment analysis
#echo -e "universe"
#col=$(head -1 ${prefix}_allhits.txt | tr -s '\t' '\n' | nl -nln |  grep "gene_id" | cut -f1 | tail -1)
#cut -f${col} ${prefix}_allhits.txt | grep -v "NA" | uniq > $universe

echo -e "Finished peak annotation:" 2>&1 | tee -a $log
