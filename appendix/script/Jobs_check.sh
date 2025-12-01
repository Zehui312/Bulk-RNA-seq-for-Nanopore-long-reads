# check_jobs_shell="basecalling.sh"
while getopts "f:l:" opt; do
    case $opt in
        f)
            check_jobs_shell="$OPTARG"
            ;;
        l)
            log_file="$OPTARG"
            ;;
        *)
            echo "Usage: $0 [-f job_script_file] [-l log_file_pattern]"
            exit 1
            ;;
    esac
done

# Set default log_file pattern if not provided

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Checking jobs in the script file: ${check_jobs_shell}"
sleep 60
echo "Log file pattern: `ls ${log_file}*.out`"
echo $submit_jobs "jobs have been submitted, waiting for successful jobs to complete..."

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

submit_jobs=$(grep "^" -c $check_jobs_shell)

waiting_time=7200
for i in `seq 1 ${waiting_time}`;do
successful_jobs=$(grep "Successfully completed" ${log_file}*.out |wc -l)

if [ $submit_jobs -eq $successful_jobs ]; then
   echo "All jobs for " ${check_jobs_shell} " have been successfully completed."
   break
else
   echo "The ${check_jobs_shell} running Current status: submitted jobs (${submit_jobs}). Successful jobs: ${successful_jobs}. Waiting for 60 seconds..."
   sleep 60
fi
done


