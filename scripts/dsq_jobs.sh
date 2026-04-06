# Notes to run a dsq job on the HPC
# First 
module load dSQ
# The 02-build_jobs.R file will generate the name 
# save this e.g. 
job_name=job_array_make_datasets_for_03-gen_data.txt
echo $job_name

# Then run, changing parameters (in the clim_risk directory)
dsq --job-file ./jobs/$job_name --mem-per-cpu 30g -t 05:00:00 --mail-type ALL --mail-type ALL --partition day --mail-user shubhi.sharma@yale.edu --out=log/slurm-%A_%a.out

# Save job id 
jobid=64641399

# Check status 
dsqa -j $jobid

# cat job_37276702_status.tsv
dsqa -j $jobid --tsv > jobs/dsq_jobs_${jobid}_status.tsv

# Additional 
# To count number of files in a directory
find /vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/spatial_pred -maxdepth 1 -type f | wc -l

# To check disk usage of a directory
du -sh /vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/spatial_pred