# Running on the Grid

## Preparing the tar file

There needs to be `irmaData.tgz` in this directory. It is constructed by,

```bash
mkdir tar ; cd tar
mkdir lib ; mkdir fcl
cp ../../project/build/lib/* lib/
cp ../../project/src/irmaData.fcl fcl/
tar cvzf ../irmaData.tgz lib/ fcl/
cd ..
```

## Preparing to run

Be sure the output file is ready (it is set in `runGridJob.sh`).

## Checking on the jobs

```bash
jobsub_q --user <USER>
```
