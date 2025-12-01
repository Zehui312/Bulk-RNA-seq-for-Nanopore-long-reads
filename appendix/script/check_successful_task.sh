grep Exited *.out |grep Job |wc -l |awk '{print "The number of failed jobs: "$1}'
grep "Successfully completed" *.out |wc -l |awk '{print "The number of successful jobs: "$1}'
ls *.err |wc -l |awk '{print "The number of all files: "$1}'
