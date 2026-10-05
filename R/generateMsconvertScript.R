#' Generates slurm script for submitting msconvert to HPC
#'
#' This function generates a slurm script to run msconvert on HPC cluster.
#' It is a helper function for the processDiannMSdap.R function.
#'
#' @param rawFiles Path where your raw files are located.
#' @param inputPath Path where input files are stored for processDiannMSdap.R function. Uses parent directory to this to save logs on a log folder.
#' @param jobName Name for the job submission. Defaults to "msconvert".
#' @param cpus Integer for how many cpus should be used for the submission. Defaults to 15.
#' @param mem How much memory should be allocated for the submission. This should be an integer followed by G, like 50G. Defaults to 100G.
#' @param time How much time should be requested for the submission. This should be integers with the following pattern: Days-Hours:Minutes:Seconds. Defaults to 1 day (01-00:00:00).
#' @return Creates slurm scripts to run MSdap and post processing on HPC.
#' @keywords internal
#' @noRd

generateMsconvertScript <- function(rawFiles,
                                  inputPath,
                                  jobName = "msconvert",
                                  cpus = 15, mem = "100G",
                                  time = "01-00:00:00") {

  log_dir <- file.path(dirname(inputPath), "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  glue::glue(r'(
#!/bin/bash
#SBATCH --job-name=<<jobName>>
#SBATCH --cpus-per-task=<<cpus>>
#SBATCH --mem=<<mem>>
#SBATCH --time=<<time>>
#SBATCH --output="<<log_dir>>/%j_%A_%a.out"
#SBATCH --error="<<log_dir>>/%j_%A_%a.err"
#SBATCH --mail-type=ALL

cd <<dirname(inputPath)>>
mkdir -p <<log_dir>>

trap "cp <<log_dir>>/$SLURM_JOB_ID* <<log_dir>>/" EXIT

# Set variables for the configs
MSCONVERT="/opt/apps/dev/containers/msconvert/1.0/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif"

# List all raw files to be converted
# Directory to search in (current directory by default)
DIR_RAW=${1:-"<<rawFiles>>"}

# Singularity msconvert use
echo "Starting msconvert..."

find "$DIR_RAW" \( -name '*.raw' -o -name '*.RAW' -o -name '*.wiff' -o -name '*.wiff2' -o -type d -name '*.d' \) | while read -r f; do

  # Check if the file is a supported file format or a .d directory
  if [[ "$f" == *.raw || "$f" == *.RAW || "$f" == *.wiff || "$f" == *.wiff2 ]]; then
    # Change the extension to .mzML
    f2="${f%.*}.mzML"

    # Skip conversion if the .mzML file already exists
    if [ ! -f "$f2" ]; then
      outdir=$(dirname "$f2")

      # Ensure the output directory exists
      mkdir -p "$outdir"

      # Run msconvert with Singularity and Wine
      singularity exec -B /home/gd2417/mywineprefix:/mywineprefix /opt/apps/dev/containers/msconvert/1.0/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif \
      mywine msconvert --64 --zlib --filter "peakPicking" --filter "zeroSamples removeExtra 1-" --outdir "$outdir" "$f"

      echo "Converted $f to $f2"
    else
      echo "Skipping $f as $f2 already exists"
    fi

  elif [[ "$f" == *.d ]]; then
    # For .d directories, just append .mzML
    f2="${f}.mzML"

    # Skip conversion if the .mzML file already exists
    if [ ! -f "$f2" ]; then
      outdir=$(dirname "$f2")

      # Ensure the output directory exists
      mkdir -p "$outdir"

      # Run msconvert with Singularity and Wine
      singularity exec -B /home/gd2417/mywineprefix:/mywineprefix /opt/apps/dev/containers/msconvert/1.0/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif \
      mywine msconvert --64 --zlib --combineIonMobilitySpectra --filter "msLevel 1-1" --outdir "$outdir" "$f"

      echo "Converted $f to $f2"
    else
      echo "Skipping $f as $f2 already exists"
    fi
  fi
done

DEST="<<inputPath>>"
mkdir -p "$DEST"
find "$DIR_RAW" -name '*.mzML' -exec cp {} "$DEST/" \;

echo "========================================="
echo "Done!"
echo "End time: $(date)"
echo "========================================="
echo "Resource usage:"
sacct -j $SLURM_JOB_ID \
--format=JobID,Elapsed,CPUTime,MaxRSS,State \
--units=G

command -v job_stats && job_stats $SLURM_JOB_ID | sed "s/\x1B\[[0-9;]*m//g"
echo "========================================="
  )', .open = "<<", .close = ">>")
}
