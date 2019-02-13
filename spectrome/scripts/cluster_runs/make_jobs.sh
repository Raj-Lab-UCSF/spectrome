num_dat=`cat MEG_filenames.dat | wc -l`

echo $num_dat

i=1
while [ $i -le $num_dat ]
do

filename=`head -n $i MEG_filenames.dat|tail -n 1`

cat job_template.sh | sed "s/NNN/${filename}/" > all_jobs/job_$i.sh

i=$(($i+1))
done
