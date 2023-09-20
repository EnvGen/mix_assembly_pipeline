#!/bin/bash -l
myfile=$2"/"$1

 mkdir -p $3/taxonomy/ && cd "$_"
 buildNCBITax=$(cat << 'EOF'
 BEGIN{ids["root"]=1;
 rank["c"]="class"; rank["d"]="superkingdom"; rank["f"]="family"; rank["g"]="genus"; rank["o"]="order"; rank["p"]="phylum";rank["s"]="species"; taxCnt=1;print "1\t|\t1\t|\tno rank\t|\t-\t|" > "nodes.dmp";
 print "1\t|\troot\t|\t-\t|\tscientific name\t|" > "names.dmp";
 }
 /^>/{ str=$2
 for(i=3; i<=NF; i++){ str=str" "$i}
 n=split(str, a, ";"); prevTaxon=1;
 for(i = 1; i<=n; i++){
 if(a[i] in ids){
 prevTaxon=ids[a[i]];
 }else{ taxCnt++;  split(a[i],b,"_"); printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", taxCnt, prevTaxon, rank[b[1]]) > "nodes.dmp";
 printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n", taxCnt, b[3]) >"names.dmp"; ids[a[i]]=taxCnt; prevTaxon=ids[a[i]]; }
 }
 gsub(">", "", $1);
 printf("%s\t%s\n", $1, ids[a[n]]) > "mapping";
 }
EOF
)
awk -F'\\[loc' '{ print $1 }' "$myfile" | awk "$buildNCBITax"
 touch merged.dmp
 touch delnodes.dmp
 cd $2
