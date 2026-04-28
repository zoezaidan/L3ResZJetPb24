#!/usr/bin/env python3
"""
HTCondor submission script for residualanalysis batch jobs

Usage:
    python3 batch/submit_condor.py \
        --era PHOTONHP_FULL \
        --input /path/to/input/directory \
        --input-type directory \
        --output-dir /path/to/output \
        --output-tag photonjet_v1 \
        --files-per-job 50 \
        --analysis PhotonJet \
        --dry-run
"""

import os
import sys
import argparse
import subprocess
from datetime import datetime

def count_files_recursive(path):
    """Count ROOT files in directory recursively"""
    count = 0
    for root, dirs, files in os.walk(path):
        count += len([f for f in files if f.endswith('.root')])
    return count

def count_files_filelist(path):
    """Count files in a filelist"""
    count = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                count += 1
    return count

def ensure_x509_proxy(proxy_path):
    """Ensure an X509 proxy exists. Create one if missing."""
    if os.path.isfile(proxy_path):
        print(f"X509 proxy: {proxy_path}")
        return proxy_path

    print(f"WARNING: X509 proxy not found at {proxy_path}")
    print("Attempting to create proxy with: voms-proxy-init --rfc --voms cms")

    cmd = ['voms-proxy-init', '--rfc', '--voms', 'cms', '--out', proxy_path]
    try:
        result = subprocess.run(cmd, text=True)
    except FileNotFoundError:
        print("Error: voms-proxy-init command not found in PATH.")
        print("Please set up the grid environment and try again.")
        sys.exit(1)

    if result.returncode != 0 or not os.path.isfile(proxy_path):
        print("Error: Failed to create X509 proxy.")
        print("Run manually: voms-proxy-init --rfc --voms cms")
        sys.exit(1)

    print(f"Created X509 proxy: {proxy_path}")
    return proxy_path

def main():
    parser = argparse.ArgumentParser(description='Submit residualanalysis batch jobs to HTCondor')
    parser.add_argument('--era', required=True, help='Era/dataset name for output naming')
    parser.add_argument('--input', required=True, help='Input path (directory/file/filelist)')
    parser.add_argument('--input-type', choices=['directory', 'file', 'filelist'],
                        default='directory', help='Input type')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--output-tag', default='batch', help='Output tag')
    parser.add_argument('--files-per-job', type=int, default=10, help='Files per job')
    parser.add_argument('--events-per-job', type=int, default=-1, help='Max events per job (-1=all)')
    parser.add_argument('--max-jobs', type=int, default=-1, help='Max number of jobs to submit (-1=all)')
    parser.add_argument('--mc', action='store_true', help='Is MC sample')
    parser.add_argument('--jet-id', action='store_true', help='Apply jet ID cuts')
    parser.add_argument('--analysis', choices=['PhotonJet', 'Dijet', 'JER'],
                        default='PhotonJet', help='Analysis type')
    parser.add_argument('--flavour', default='workday',
                        choices=['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow'],
                        help='HTCondor job flavour')
    parser.add_argument('--jet-tree', default='ak4PFJetAnalyzerSDZcut1/t',
                        help='Jet tree path (e.g. ak4PFJetAnalyzer/t, ak4PFJetAnalyzerSDZcut1/t)')
    parser.add_argument('--dry-run', action='store_true', help='Create files but do not submit')

    args = parser.parse_args()

    # Setup paths
    base_dir = '/eos/home-b/bharikri/lxplus_private/EGamma/residualanalysis'
    analysis_dir = os.path.join(base_dir, 'fillhistograms')
    cmssw_base = '/eos/home-b/bharikri/lxplus_private/EGamma/CMSSW_15_1_0_patch3/src'

    # Find or create X509 grid proxy for remote file access
    x509_proxy = os.environ.get('X509_USER_PROXY', f'/tmp/x509up_u{os.getuid()}')
    x509_proxy = ensure_x509_proxy(x509_proxy)

    # Count total files
    if args.input_type == 'directory':
        total_files = count_files_recursive(args.input)
    elif args.input_type == 'filelist':
        total_files = count_files_filelist(args.input)
    else:
        total_files = 1

    # Calculate number of jobs
    total_jobs = (total_files + args.files_per_job - 1) // args.files_per_job
    
    # Apply max jobs limit if specified
    if args.max_jobs > 0 and total_jobs > args.max_jobs:
        total_jobs = args.max_jobs
        print(f"WARNING: Limiting to {total_jobs} jobs (use --max-jobs to change)")

    print(f"Input: {args.input}")
    print(f"Total files: {total_files}")
    print(f"Files per job: {args.files_per_job}")
    print(f"Total jobs: {total_jobs}")

    # Create batch directory with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    batch_dir = os.path.join(base_dir, 'batch/jobs/', f'{args.era}_{timestamp}')
    os.makedirs(batch_dir, exist_ok=True)
    os.makedirs(os.path.join(batch_dir, 'logs'), exist_ok=True)

    print(f"Batch directory: {batch_dir}")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Create job wrapper script
    wrapper_path = os.path.join(batch_dir, 'run_job.sh')
    wrapper_content = f'''#!/bin/bash
# Job wrapper for residualanalysis batch

ANALYSIS=$1
INPUT_PATH=$2
OUTPUT_TAG=$3
IS_MC=$4
JET_ID=$5
INPUT_TYPE=$6
MAX_FILES=$7
MAX_EVENTS=$8
OUTPUT_DIR=$9
BATCH_INDEX=${{10}}
TOTAL_BATCHES=${{11}}
JET_TREE=${{12}}

echo "Starting job at $(date)"
echo "Analysis: $ANALYSIS"
echo "Input: $INPUT_PATH"
echo "Jet ID: $JET_ID"
echo "Jet tree: $JET_TREE"
echo "Batch: $BATCH_INDEX of $TOTAL_BATCHES"

# Setup X509 proxy for xrootd access
export X509_USER_PROXY=${{X509_USER_PROXY:-{x509_proxy}}}
echo "X509 proxy: $X509_USER_PROXY"

# Setup CMSSW
cd {cmssw_base}
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

# Go to analysis directory
cd {analysis_dir}

# Run analysis
root -l -b -q "analyse_${{ANALYSIS}}.cc(\\"$INPUT_PATH\\", \\"$OUTPUT_TAG\\", $IS_MC, $JET_ID, \\"$INPUT_TYPE\\", $MAX_FILES, $MAX_EVENTS, \\"$OUTPUT_DIR\\", $BATCH_INDEX, $TOTAL_BATCHES, \\"$JET_TREE\\")"

echo "Job completed at $(date)"
'''

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)
    os.chmod(wrapper_path, 0o755)

    # Create job arguments file
    args_file = os.path.join(batch_dir, 'job_args.txt')
    is_mc = 'true' if args.mc else 'false'
    jet_id = 'true' if args.jet_id else 'false'

    with open(args_file, 'w') as f:
        for i in range(total_jobs):
            line = f"{args.analysis} {args.input} {args.output_tag} {is_mc} {jet_id} {args.input_type} {args.files_per_job} {args.events_per_job} {args.output_dir} {i} {total_jobs} {args.jet_tree}\n"
            f.write(line)

    # Create condor submit file
    submit_path = os.path.join(batch_dir, 'condor.sub')
    submit_content = f'''# HTCondor submit file for residualanalysis
# Generated: {datetime.now().isoformat()}

universe = vanilla
executable = {wrapper_path}
arguments = $(args)

output = {batch_dir}/logs/job_$(Process).out
error = {batch_dir}/logs/job_$(Process).err
log = {batch_dir}/logs/condor.log

should_transfer_files = NO
x509userproxy = {x509_proxy}
use_x509userproxy = true
+JobFlavour = "{args.flavour}"
request_cpus = 1
request_memory = 4000

queue args from {args_file}
'''

    with open(submit_path, 'w') as f:
        f.write(submit_content)

    print(f"\nCreated {total_jobs} job configurations")
    print(f"Submit file: {submit_path}")

    if not args.dry_run:
        print("\nSubmitting to HTCondor...")
        result = subprocess.run(['condor_submit', submit_path],
                               capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            sys.exit(1)
        print("\nJobs submitted successfully!")
        print(f"Monitor with: condor_q")
        print(f"Logs in: {batch_dir}/logs/")
    else:
        print("\n[DRY RUN] Would submit:")
        print(f"  condor_submit {submit_path}")
        print("\nTo submit manually:")
        print(f"  condor_submit {submit_path}")

if __name__ == '__main__':
    main()
