for folder in $(find . -maxdepth 1 -mindepth 1 -type d); do
    # Get the base name. 
    base_name=$(basename "$folder")
    # echo "${base_name}"

    # Create (command to create) the zip file. 
    # echo "zip -r ${base_name}.zip ${base_name}/"
    zip -r ${base_name}.zip ${base_name}/
done # > zip_commands 
# n_jobs=$(find . -maxdepth 1 -mindepth 1 -type d | wc -l)
# split -n l/${n_jobs} --numeric-suffixes zip_commands job.zip.folder.
# for i in `ls job.zip.folder.*`; do cat $i | bsub -o $i.lsf.out -e $i.lsf.err -q e80medium -W 6:00 -M 80 -R 'rusage[mem=80]' -n 28 -u lwfong@mdanderson.org -J $i -P zip_folders; done