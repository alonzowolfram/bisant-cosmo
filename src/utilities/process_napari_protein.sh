#!/bin/bash

# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/napari-protein-stitching/napari-cosmx-protein-stitching.html

# Flags
# -i: Path containing folders, each folder = 1 sample. -i, -f for stitch-images will be derived from this.
# -o: Path to which to save zarr output. 
# https://www.baeldung.com/linux/use-command-line-arguments-in-bash-script
while getopts i:o:n: flag
do
    case "${flag}" in
        i) input_folder=${OPTARG};;
        o) output_folder=${OPTARG};;
        n) run_name=${OPTARG};;
    esac
done 
# echo "Input folder: $input_folder"
# echo "Output folder: $output_folder"

# Loop for each folder given by -i.
for folder in ${input_folder}/*; do
# $folder will be the complete path. 
    
    # Check if it's a directory.
    if [ -d "$folder" ]; then

        # Extract the basename.
        base_name=$(basename "$folder")
        # Create the output folder.
        individual_output_folder=${output_folder}/${base_name}
        mkdir ${individual_output_folder}

        # Create commands to run stitch-images and stitch-expression.
        # First we need to get the folder that's not the Logs folder inside $folder.
        raw_dat_folder=$(find -L $folder -maxdepth 1 -mindepth 1 -type d | grep -v "Logs") # Previously piped grep into `sed "2q;d"``
        # echo "Raw data folder: ${raw_dat_folder}"
        # $raw_dat_folder will be the complete path. 
        # Next we need to get the folder inside AnalysisResults.
        # This assumes we have only one folder inside AnalysisResults,
        # and takes the first folder returned by the expression `find $folder -type d | sed '2q;d'`
        voting_folder=$(find -L $raw_dat_folder/AnalysisResults/ -type d | sed '2q;d')
        # echo "Voting folder: ${voting_folder}"

        # Echo the commands to run the programs.
        # echo "CellStatsDir folder: ${raw_dat_folder}/CellStatsDir/"
        # echo "RunSummary folder: ${raw_dat_folder}/RunSummary/"

        echo "source /rsrch6/home/genomic_med/lwfong/miniconda3/etc/profile.d/conda.sh; conda activate napari-cosmx; export PATH=~/miniconda3/envs/napari-cosmx/bin:$PATH; which stitch-images; which stitch-expression; stitch-images -i ${raw_dat_folder}/CellStatsDir/ -f ${raw_dat_folder}/RunSummary/ -o ${individual_output_folder}; stitch-expression -i ${raw_dat_folder} -o ${individual_output_folder}; exit_code=\$?; if [ \$exit_code -eq 0 ]; then echo \"Napari job ${run_name}, folder ${folder}, individual output folder ${individual_output_folder} completed with exit code 0. Results available at ${individual_output_folder}.\" | mail -s \"Napari job ${run_name}, folder ${folder}, individual output folder ${individual_output_folder} completed successfully\" lwfong@mdanderson.org; else echo \"Napari job ${run_name}, folder ${folder}, individual output folder ${individual_output_folder} completed with non-zero exit code \${exit_code}\" | mail -s \"Napari job ${run_name}, folder ${folder}, individual output folder ${individual_output_folder} failed\" lwfong@mdanderson.org; fi"
    fi
done > Commands # Write to Commands file.
# Split Commands file into 1 job each. 
n_jobs=$(find -L $input_folder -maxdepth 1 -mindepth 1 -type d | wc -l) # https://stackoverflow.com/questions/11071442/bash-find-exclude-parent
# Create a random string for this job so it doesn't overlap with other jobs.
random_string=$(echo $RANDOM$RANDOM | sha256sum | head -c 16)
echo "Total lines in Commands: $(wc -l < Commands)"
echo "Number of jobs: $n_jobs"
split -l $(( $(wc -l < Commands) / n_jobs )) --numeric-suffixes Commands job.process.napari.${random_string}. # split -n l/${n_jobs} --numeric-suffixes Commands job.process.napari.${random_string}. # l is needed to prevent from splitting a line in two, per the manual.
for i in `ls job.process.napari.${random_string}*`; do cat $i | bsub -o ${output_folder}/$i.lsf.out -e ${output_folder}/$i.lsf.err -q medium -W 4:00 -M 80 -R 'rusage[mem=80]' -n 28 -u lwfong@mdanderson.org -J $i -P process_napari; done

# bsub -Is -q medium -W 4:00 -M 80 -R 'rusage[mem=80]' -n 28 /usr/bin/bash