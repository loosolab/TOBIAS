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
	echo -e "	-i <string> []			input bed file for annotation"
	echo -e "	-o <string> []			output bed file with path"
	echo -e "	-g <string> []			gtf file for annotation, make sure filter attribute suits"

	echo -e "Optional"
	echo -e "	-p <string> [basename bed _header.txt]	header file"
	echo -e "	-l <string> [basename bed .log]	log file"
	echo -e "	-c <int> [1]			number of cores that should be used"
	echo -e "	-t <string> [gene]		feature for annotation"
	echo -e "	-r <string> [start]		feature anchor for distance calculation"
	echo -e "	-e <string> [[10000,1000]]	distance to feature, either one distance or distance interval"
	echo -e "	-a <string> [gene_biotype]	attribute for filtering, attention 'gene_biotype' in ensemble, but 'gene_type' in gencode"
	echo -e "	-w <string> [protein_coding]	attribute value for filtering"
	echo -e "	-n <string> [gene_name,gene_id, ]	show attributes for uropa output"
	echo -e "	-v				print stdout only to log file"
	exit 1
}
export LC_ALL=C

# defaults
cores=1
feature="gene"
feature_anchor="start"
distance=[10000,1000]
#filter_attribute="gene_biotype"
#value_attribute="protein_coding"
show_attribute='"gene_biotype", "gene_name", "gene_id"'

# no option given
if [ $# -eq 0 ]; then
	usage
	exit -1
fi

# Process command line options
progname=${0##*/}
shortopts="i:o:j:p:l:c:g:t:r:e:a:w:n:vh"
longopts="bed:,annotation:,universe:,header:,log:,cores:,gtf:,feature:,feature_anchor:,distance:,filter_attribute:,attribute_value:,show_attribute:,verbose"
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
        ( -v|--verbose          ) verbose=yes ;;
    esac
    shift
done

echo -e "feature: $feature \ndistance: $distance \nshow_attribute: $show_attribute"

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

if [[ ! -z $filter_attribute ]] && [[ ! -z $value_attribute ]]; then
	cat > $json << EOL
{
	"queries":[
	{"feature":"$feature", "feature.anchor":"$feature_anchor", "distance":$distance, "filter.attribute":"$filter_attribute", "attribute.value":"$value_attribute", "show.attributes":[$show_attribute]}
	],
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
# tr '**' '\t' |
echo -e "header"
#echo -e "tf_chr\ntf_start\ntf_end\ntf_name\ntf_score\ntf_strand" > $header
attributes=$(($(head -1 ${prefix}_finalhits.txt | awk '{print NF}')-1))
`head -1 ${prefix}_finalhits.txt | awk 'BEGIN { OFS = "\n" } ; { print $2,$3,$5,$1,$6,$7,$8,$9,$10,$11 }' - >> $header`
`head -1 ${prefix}_finalhits.txt | cut -f15-${attributes} - | sed -e 's/\t/\n/g' >> $header`
#echo -e "fp_score" >> $header
`awk 'BEGIN { OFS = "\t" } ; { print $2,$3,$5,$1,$6,$7,$8,$9,$10,$11 }' ${prefix}_finalhits.txt | sed -e 1d |  sed -e 's/\t\t/\t/g' > $annotation'.tmp'` 	
cut -f15-${attributes} ${prefix}_finalhits.txt | sed -e 1d | paste -d"\t" $annotation'.tmp' - > $annotation
sort -k1,1 -k2,2n $annotation > $annotation'.tmp'
mv $annotation'.tmp' $annotation
rm $annotation'.tmp'


echo -e "universe"
col=$(head -1 ${prefix}_allhits.txt | tr -s '\t' '\n' | nl -nln |  grep "gene_id" | cut -f1 | tail -1)
cut -f${col} ${prefix}_allhits.txt | grep -v "NA" | uniq > $universe

echo -e "Finished peak annotation:" 2>&1 | tee -a $log
