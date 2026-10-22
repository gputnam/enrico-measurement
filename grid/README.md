# Grid production of makedf dataframes

Requirements: `jobsub_lite` with a valid token (`htgettoken -a htvaultprod.fnal.gov -i sbnd`),
membership in the `sbnd` group.  Jobs run in the `fnal-wn-el9` container with spack ROOT from
cvmfs, stream the MINERvA Open Data over xrootd and write to
`/pnfs/sbnd/scratch/users/$USER/MvA-Loop/df/<tag>/<sample>/`.

| script | purpose |
|---|---|
| `build_tarball.sh <tag>` | build makedf, package binary + column configs + aux files |
| `make_joblists.py` | one MC file per job; data files grouped to ~5 GB per job |
| `submit.sh <tag> <joblist>` | jobsub_submit (-N = number of lines) |
| `gridjob.sh` | worker: runs line `$PROCESS` of the job list, copies output + log with ifdh |
| `production.sh <tag> [samples]` | all of the above for a full production |
| `check_outputs.py <joblist> <outdir>` | completeness / POT check, writes a resubmission list |

Typical full production:
```
grid/production.sh v1 me_fhc_mc me_fhc_data
grid/check_outputs.py joblists/v1/me_fhc_mc.txt /pnfs/sbnd/scratch/users/$USER/MvA-Loop/df/v1/me_fhc_mc --resubmit-list joblists/v1/me_fhc_mc_resub.txt
OUTNAME=me_fhc_mc grid/submit.sh v1 joblists/v1/me_fhc_mc_resub.txt   # same output dir
```
A reprocessing with new columns but the same event sample reuses the previous tag's job lists
(v2 used `joblists/v1/*.txt` unchanged: `mkdir -p joblists/v2 && cp joblists/v1/me_fhc_{mc,data}.txt joblists/v2/`,
then `grid/build_tarball.sh v2` and `grid/submit.sh v2 joblists/v2/<sample>.txt`), so POT and file
lists are identical; compare with `grid/check_outputs.py` on both tags.

Throughput: ~80 s per 21 GB MC file and ~30 s per data file on fast dCache pools (some pools are
much slower); output ~70 MB (v1) / ~80 MB (v2) per MC file, ~10 / ~12 MB per data job.

Inputs on an unavailable dCache pool make makedf exit with rc 6 after processing the other files.
For such persistent failures resubmit with `-e MAKEDF_ALLOW_MISSING=1` (extra jobsub argument) to keep
the partial output; POT bookkeeping stays correct because `meta` lists only processed files, and
the missing inputs are listed in the job log.
