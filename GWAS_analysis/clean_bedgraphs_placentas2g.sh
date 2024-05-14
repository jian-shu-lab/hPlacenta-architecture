bedops_cell=/data/deyk/extras/BEDOPS/bin
bedtools_cell=/data/deyk/extras/bedtools2/bin
bed_cell=~/data/BEDFILES/SHAREseq_S2G


TASKFILE=~/data/placentas2g.txt

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   annot_name=`echo $line | awk '{print $1}'`
   input_cell=$bed_cell/$annot_name
   echo  $input_cell
   names=`ls $input_cell`
   for name in $names
   do
       $bedtools_cell/bedtools sort -i $input_cell/$name > $input_cell/$name.2
       $bedtools_cell/bedtools merge -i $input_cell/$name.2 -c 4 -o max > $input_cell/$name.3
       mv $input_cell/$name.3 $input_cell/$name
       rm $input_cell/$name.2
   done
done 
