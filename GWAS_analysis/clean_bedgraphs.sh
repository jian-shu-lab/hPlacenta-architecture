bedops_cell=/data/deyk/extras/BEDOPS/bin
bedtools_cell=/data/deyk/extras/bedtools2/bin
bed_cell=~/data/BEDFILES/placenta


TASKFILE=~/data/placenta.txt

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   annot_name=`echo $line | awk '{print $1}'`
   input_cell=$bed_cell/$annot_name
   echo  $input_cell
   names=`ls $input_cell | cut -f 1 -d '.'`
   for name in $names
   do
       $bedtools_cell/bedtools sort -i $input_cell/$name.bed > $input_cell/$name.2.bed
       $bedtools_cell/bedtools merge -i $input_cell/$name.2.bed -c 4 -o max > $input_cell/$name.3.bed
       mv $input_cell/$name.3.bed $input_cell/$name.bed
       rm $input_cell/$name.2.bed
   done
done 
