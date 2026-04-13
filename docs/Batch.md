# Batch

Batch mode is only for histogram production. Merge first, then run the L2, JER, or L3 derivation and fit macros on the merged ROOT file.

## Minimal flow

| Step | Helper | Result |
| --- | --- | --- |
| submit jobs | `batch/submit_condor.py` | many histogram-filler jobs |
| merge shards | `batch/merge_outputs.sh` | one merged ROOT file |
| validate locally | `batch/test_batch.sh` | dry-run submission check |

## Submission Interface

`batch/submit_condor.py` is the user-facing entry point. Its command-line arguments are:

| Argument | Meaning |
| --- | --- |
| `--era` | dataset or era label used in naming |
| `--input` | directory, file, or filelist path |
| `--input-type` | one of `directory`, `file`, or `filelist` |
| `--output-dir` | directory where shard outputs are written |
| `--output-tag` | output tag appended to shard filenames |
| `--files-per-job` | number of input files per Condor job |
| `--events-per-job` | maximum events per job, `-1` means all |
| `--max-jobs` | cap on the number of submitted jobs, `-1` means all |
| `--mc` | mark the sample as MC |
| `--jet-id` | enable jet ID in the histogram macro |
| `--analysis` | one of `PhotonJet`, `Dijet`, or `JER` |
| `--flavour` | HTCondor runtime flavour |
| `--jet-tree` | jet tree path such as `ak4PFJetAnalyzer/t` or `ak4PFJetAnalyzerSDZcut1/t` |
| `--dry-run` | create job files but do not submit |

## Typical submission

```bash
python3 batch/submit_condor.py \
  --era PHOTONHP \
  --input batch/input/filelist_forests_2024ppRef_HP_ALL.txt \
  --input-type filelist \
  --output-dir /path/to/output \
  --output-tag photonjet \
  --files-per-job 50 \
  --analysis PhotonJet
```

For the dijet and JER branches, only the `--analysis`, `--mc`, `--jet-id`, and input naming change. The batch layer still only produces histogram ROOT files.

## Merge

```bash
./batch/merge_outputs.sh /path/to/output PHOTONHP photonjet
```

The merged ROOT file from batch processing is the input to the downstream derivation macros:

| Workflow | Next macro after merge |
| --- | --- |
| L2 residuals | `L2Residual/deriveL2_from3D.C` |
| JER truth or scale factors | `JER/MCJER.C`, `JER/MCJPR.C`, `JER/JERSF_RMS.C`, or `JER/JERSF_fits.C` |
| L3 residuals | `L3Residual/deriveL3_from_photonjet.C` |

## Notes

| Item | Meaning |
| --- | --- |
| input modes | same `era`, `file`, `directory`, and `filelist` modes already supported by the analysis macros |
| filelists | stored under `batch/input/` |
| naming | keep output tags stable so MC and data files are easy to pair later |