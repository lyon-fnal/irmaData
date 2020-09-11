# Running on the Grid

## Preparing the tar file

There needs to be `irmaData.tgz` in this directory. It is constructed by,

```bash
mkdir tar ; cd tar
mkdir lib ; mkdir fcl
cp ../../project/build/lib/* lib/
cp ../../project/src/irmaData.fcl fcl/
tar cvzf ../irmaData.fcl lib/ fcl/
cd ..
```

## Checking on the jobs

```bash
jobsub -q --user <USER>
```
