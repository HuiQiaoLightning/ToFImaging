#!/bin/bash
#Read the jobnum
jobnum=$(<job.txt)

#Write the job
newjobnum=$((jobnum+1))
echo "$newjobnum" > job.txt

#Start the job
echo "Starting the job $jobnum"
matlab -nojvm -logfile "matlog_job_$jobnum.log" -r "reconstruct_video_fitmodel_launch($jobnum)"

