annot_cell=~/data/ANNOTATIONS_hg38/SHAREseq_S2G
output_cell_pre=[OUTPUT_PATH]/placenta/
baseline_version=[BASELINE_VERSION]
output_cell=$output_cell_pre/$baseline_version

sumstats_taskfile=~/data/traits_placenta.txt
#sumstats_taskfile=/n/groups/price/kushal/singlecellLDSC/data/traits_bio.txt

IFS="
"

module load gcc/10.2.0
module load R

flag=0
index_in_results=1 ## which annotation to choose from the .results file in case of multiple annotations


for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
do
sumstats_file=`echo $step | awk '{print $1}'`
echo $sumstats_file
sumstats_file2=${sumstats_file%.sumstats}

counter1=0
for step2 in `ls $output_cell | awk '{print $1}' | sort | uniq`;
do
    annot_name=`echo $step2 | awk '{print $1}'`
    if [ ! -f $output_cell/$annot_name/${sumstats_file2}_ldsc_postprocess.txt ]
    then
	counter1=$(($counter1+1))
    fi
done

if (( $counter1 > 0 ))
then
    echo $sumstats_file2
    cmd="Rscript ldsc_postprocess.R  $annot_cell $output_cell $sumstats_file $flag $index_in_results"
    bsub -W 270 -R "rusage[mem=20]" -e ldsc_post.err -o ldsc_post.out -n 1 "$cmd"
   # sbatch --time=40:00 --mem=20000 --output=ldsc_post.out --error=ldsc_post.err -p short -c 1 --wrap="$cmd"
fi
done
