# Job submission script

jobsub_submit -N 300 -G gm2 -e ERA=2C --dataset_definition=gm2pro_daq_offline_dqc_20p_run2C --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh

sleep 30

jobsub_submit -N 250 -G gm2 -e ERA=2D --dataset_definition=gm2pro_daq_offline_dqc_run2D --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh

sleep 30

jobsub_submit -N 200 -G gm2 -e ERA=2E --dataset_definition=gm2pro_daq_offline_dqc_run2E --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh

sleep 30

jobsub_submit -N 200 -G gm2 -e ERA=2F --dataset_definition=gm2pro_daq_offline_dqc_run2F --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh

sleep 30

jobsub_submit -N 50 -G gm2 -e ERA=2G --dataset_definition=gm2pro_daq_offline_dqc_run2G --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh

sleep 30

jobsub_submit -N 80 -G gm2 -e ERA=2H --dataset_definition=gm2pro_daq_offline_dqc_run2H --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=15h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh
