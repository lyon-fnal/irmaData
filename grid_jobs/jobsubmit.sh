# Job submission script

jobsub_submit -N 250 -G gm2 --dataset_definition=gm2pro_daq_offline_dqc_run2E --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --cpu=1 --memory=1990MB --disk=20GB --expected-lifetime=10h --role=Analysis --tar_file_name dropbox://$PWD/irmaData.tgz --use-cvmfs-dropbox file:///gm2/app/users/lyon/sl7/irmaData/grid_jobs/runGridJob.sh
