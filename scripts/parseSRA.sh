
"""
# Determine The intersection of sequencing types for humans and mice
comm -12 <(cut -f9,10 SRAdata.txt | sort | uniq -c | grep "Homo_sapien" | awk '{if($1 > 1000) {print $0}}' | cut -f2 | sort) <(cut -f9,10 SRAdata.txt | sort |
	uniq -c | grep "Mus_musculus" | awk '{if($1 > 1000) {print $0}}' | cut -f2 | sort);
"""

Species=(
    "Homo_sapiens"
    "Mus_musculus"
	"Caenorhabditis_elegans"
	"Rattus_norvegicus"
	"Saccharomyces_cerevisiae"
	"Arabidopsis_thaliana"
	"Drosophila_melanogaster"
)
Library=(
	"ChIP-Seq"
	"RNA-Seq"
	"WGS"
)

for Spec in ${Species[@]}; do
echo $Spec ;
cat SRAdata.txt | cut -f1,9,10,19 | grep "public" | grep -w "$Spec" - > species.txt

for Lib in ${Library[@]}; do
echo $Lib ;
cat species.txt | grep -w "$Lib" - | cut -f1 | shuf | head -1000 > library.txt
cat library.txt | awk -v spec=$Spec -v lib=$Lib 'OFS="," { print spec,lib,$0 }' >> all-accessions.txt;
done;
done;
