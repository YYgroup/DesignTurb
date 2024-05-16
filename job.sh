echo "#!/bin/sh">submit.sh
echo "yhrun -N 6 -n 128 -p work_short_job ./tube >out 2>&1">>submit.sh
chmod 700 submit.sh 
yhbatch -N 6 -n 128 -p work_short_job -J turb.swy  ./submit.sh
rm submit.sh
