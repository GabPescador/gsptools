#' Generates slurm script for submitting msconvert to HPC
#'
#' This function generates a slurm script to convert raw vendor files (.raw,
#' .wiff/.wiff2, .d) to .mzML via msconvert (Singularity + Wine), then copies
#' the resulting .mzML files into the MSDAP input directory.
#' It is a helper function for the processDiannMSdap.R function and MUST run
#' BEFORE generateMsdapScript.R, since MSDAP's QC report expects the .mzML
#' files to sit alongside report.tsv in `inputPath`.
#'
#' @param rawFiles Path where your raw files are located.
#' @param inputPath Path where input files are stored for processDiannMSdap.R function. Uses parent directory to this to save logs on a log folder.
#' @param jobName Name for the job submission. Defaults to "msconvert".
#' @param cpus Integer for how many cpus should be used for the submission. Defaults to 15.
#' @param mem How much memory should be allocated for the submission. This should be an integer followed by G, like 50G. Defaults to 100G.
#' @param time How much time should be requested for the submission. This should be integers with the following pattern: Days-Hours:Minutes:Seconds. Defaults to 1 day (01-00:00:00).
#' @return Creates slurm script to convert raw files and stage .mzML for MSdap.
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

# --- Config ---
MSCONVERT="/opt/apps/dev/containers/msconvert/1.0/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif"
DIR_RAW=${1:-"<<rawFiles>>"}
DIR_MZML=<<inputPath>>

# --- Pre-check: skip the entire find/convert/stage block if every raw file
#     already has a staged mzML counterpart in 4_input. Avoids re-running
#     msconvert on a job that's already fully done. ---
RAW_COUNT=$(find "$DIR_RAW" \( -name '*.raw' -o -name '*.RAW' -o -name '*.wiff' -o -name '*.wiff2' -o -type d -name '*.d' \) | wc -l)
STAGED_COUNT=$(find "$DIR_MZML" -maxdepth 1 -name '*.mzML' | wc -l)
 
if [ "$RAW_COUNT" -gt 0 ] && [ "$STAGED_COUNT" -ge "$RAW_COUNT" ]; then
  echo "Found $RAW_COUNT raw files and $STAGED_COUNT mzML already staged in $DIR_MZML - skipping conversion entirely."
else
  echo "Starting msconvert ($RAW_COUNT raw files found, $STAGED_COUNT already staged)..."

# --- Convert every supported raw format found under DIR_RAW ---
echo "Starting msconvert..."

# Manifest of exactly the mzML paths this run's raw/wiff/.d files map to -
# used below to stage only those files, so any unrelated mzML sitting under
# DIR_RAW (e.g. *_diatracer.mzML from another search engine) never gets
# copied into the MSDAP input folder just because it matches *.mzML.
MANIFEST=$(mktemp)

find "$DIR_RAW" \( -name '*.raw' -o -name '*.RAW' -o -name '*.wiff' -o -name '*.wiff2' -o -type d -name '*.d' \) | while read -r f; do

  # Mirror the input's subpath under DIR_RAW into DIR_MZML, so same-named
  # files from different subfolders of DIR_RAW don't collide/overwrite.
  relpath="${f#$DIR_RAW}"

  if [[ "$f" == *.raw || "$f" == *.RAW || "$f" == *.wiff || "$f" == *.wiff2 ]]; then
    # Thermo/Sciex: swap extension for .mzML, rooted in DIR_MZML
    f2="$DIR_MZML/${relpath%.*}.mzML"

    if [ ! -f "$f2" ]; then
      outdir=$(dirname "$f2")
      mkdir -p "$outdir"
      singularity exec -B /home/gd2417/mywineprefix:/mywineprefix "$MSCONVERT" \
        mywine msconvert --64 --zlib --filter "msLevel 1-1" --outdir "$outdir" "$f"
      echo "Converted $f to $f2"
    else
      echo "Skipping $f as $f2 already exists"
    fi
    echo "$f2" >> "$MANIFEST"

  elif [[ "$f" == *.d ]]; then
    # Bruker: append .mzML to the directory name, rooted in DIR_MZML
    f2="$DIR_MZML/${relpath}.mzML"

    if [ ! -f "$f2" ]; then
      outdir=$(dirname "$f2")
      mkdir -p "$outdir"
      singularity exec -B /home/gd2417/mywineprefix:/mywineprefix "$MSCONVERT" \
        mywine msconvert --64 --zlib --combineIonMobilitySpectra --filter "msLevel 1-1" --outdir "$outdir" "$f"
      echo "Converted $f to $f2"
    else
      echo "Skipping $f as $f2 already exists"
    fi
    echo "$f2" >> "$MANIFEST"
  fi
done

# --- Stage only the mzML files converted above next to report.tsv, so MSDAP
#     QC can find them (unrelated mzML files elsewhere under DIR_RAW are left
#     untouched) ---
DEST="<<inputPath>>"
mkdir -p "$DEST"
while IFS= read -r f2; do
  dest_file="$DEST/$(basename "$f2")"
  if [ -f "$f2" ]; then
    if [ -f "$dest_file" ]; then
      echo "Skipping staging of $f2 as $dest_file already exists"
    else
      cp "$f2" "$DEST/"
    fi
  fi
done < "$MANIFEST"
rm -f "$MANIFEST"
fi

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
