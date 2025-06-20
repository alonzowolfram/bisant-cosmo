#!/bin/bash

# stitch-images -i CellStatsDir/ -f RunSummary/ -o ...
# read-targets -o ... AnalysisResults/.../
# i: Path to CellLabels and morphology images. Input to stitch-images.
# f: Path to latest.fovs.csv directory. Input to stitch-images.
# o: Where to create zarr output. Input to stitch-images.
# folder: voting folder. Will be named AnalysisResults/.+ Input to read-targets.

# Flags
# -i: Path containing folders, each folder = 1 sample. -i, -f for stitch-images and the folder positional argument for read-targets will be derived from this.
# -o: Path to which to save zarr (stitch-images) and hdf5 (read-targets) output. 
# https://www.baeldung.com/linux/use-command-line-arguments-in-bash-script
while getopts i:o: flag
do
    case "${flag}" in
        i) input_folder=${OPTARG};;
        o) output_folder=${OPTARG};;
    esac
done 
# echo "Input folder: $input_folder"
# echo "Output folder: $output_folder"

python /rsrch6/home/genomic_med/lwfong/PRIME-TR/projects/CosMx-P001_bisant-cosmo/src/utilities/stitch_napari.py $input_folder $output_folder