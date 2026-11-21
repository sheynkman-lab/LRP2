#!/usr/bin/env bash
# Prepares a GitHub-hosted runner's disks for an integration test. Used by every job in
# integration_tests.yml, so the three jobs can't drift apart.
#
# The runner cannot hold the human reference, task work dirs, every container image and (for
# Apptainer) a layer cache at once ("no space left on device"). This script:
#   1. deletes preinstalled software the tests never use;
#   2. puts the Nextflow work dir and Singularity/Conda storage on whichever of / and /mnt has more
#      room, and exports the paths through $GITHUB_ENV;
#   3. optionally adds swap: pbmm2 indexing the full human genome (--preset ISOSEQ) needs more than
#      the runner's 16 GB of RAM;
#   4. for Docker, starts ci_docker_disk_guard.sh in the background. Docker keeps every image it has
#      pulled until the job ends, and the images of one full run add up to more than the disk holds.
#
# Environment:
#   CONTAINER_ENGINE  docker | singularity | conda   (default: docker)
#   CI_SWAP_GB        size of the swap file in GB, 0 for none   (default: 16)
set -euo pipefail

ENGINE="${CONTAINER_ENGINE:-docker}"
SWAP_GB="${CI_SWAP_GB:-16}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

avail_gb() { df -BG --output=avail "$1" 2> /dev/null | tail -n 1 | tr -dc '0-9'; }

echo "Free space before cleanup: / = $(avail_gb /) GB"

# Preinstalled SDKs and tools. Java in /opt/hostedtoolcache is kept (setup-java uses it), and so is
# /usr/share/miniconda. Everything here is optional: a path that doesn't exist is simply skipped.
sudo rm -rf \
  /usr/share/dotnet /usr/local/lib/android /opt/ghc /usr/local/.ghcup \
  /opt/hostedtoolcache/CodeQL /opt/hostedtoolcache/go /opt/hostedtoolcache/PyPy \
  /opt/hostedtoolcache/Ruby /opt/hostedtoolcache/node /opt/hostedtoolcache/Python \
  /usr/local/share/boost /usr/share/swift /usr/local/share/powershell /usr/local/share/chromium \
  /opt/microsoft /opt/az /opt/google /opt/pipx /usr/lib/google-cloud-sdk /usr/local/share/vcpkg \
  /usr/local/lib/node_modules /usr/local/aws-cli /usr/local/aws-sam-cli /usr/lib/firefox \
  /usr/share/kotlinc /usr/share/sbt || true
sudo rm -rf /usr/local/julia* /usr/share/az_* /usr/share/gradle-* /usr/lib/llvm-* || true
# Images, stopped containers, volumes and build cache that came with the runner
sudo docker system prune --all --volumes --force > /dev/null 2>&1 || true
sudo apt-get clean || true
sudo rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* || true

ROOT_FREE=$(avail_gb /)
MNT_FREE=0
if mountpoint -q /mnt; then MNT_FREE=$(avail_gb /mnt); fi
echo "Free space after cleanup:  / = ${ROOT_FREE} GB, /mnt = ${MNT_FREE} GB"

if (( MNT_FREE > ROOT_FREE )); then
  SCRATCH=/mnt/ci;               SCRATCH_FREE=$MNT_FREE;  OTHER=/;    OTHER_FREE=$ROOT_FREE
else
  SCRATCH="${RUNNER_TEMP}/ci";   SCRATCH_FREE=$ROOT_FREE; OTHER=/mnt; OTHER_FREE=$MNT_FREE
fi
sudo mkdir -p "${SCRATCH}/work" "${SCRATCH}/images" "${SCRATCH}/apptainer-cache" "${SCRATCH}/apptainer-tmp"
sudo chown -R "$(id -u):$(id -g)" "$SCRATCH"
{
  echo "CI_SCRATCH=${SCRATCH}"
  echo "NXF_WORK=${SCRATCH}/work"
  echo "NXF_SINGULARITY_CACHEDIR=${SCRATCH}/images"
  echo "APPTAINER_CACHEDIR=${SCRATCH}/apptainer-cache"
  echo "APPTAINER_TMPDIR=${SCRATCH}/apptainer-tmp"
  # Without this, every image is stored twice: as OCI layers in the cache and as a SIF file
  echo "APPTAINER_DISABLE_CACHE=true"
} >> "$GITHUB_ENV"
echo "Scratch dir (work dir, Singularity images, Conda envs): ${SCRATCH}"

if (( SWAP_GB > 0 )); then
  SWAP_FILE=""
  if (( OTHER_FREE >= SWAP_GB + 24 )); then
    SWAP_FILE="${OTHER%/}/ci-swapfile"
  elif (( SCRATCH_FREE >= SWAP_GB + 44 )); then
    SWAP_FILE="${SCRATCH}/ci-swapfile"
  fi
  if [[ -n "$SWAP_FILE" ]]; then
    sudo fallocate -l "${SWAP_GB}G" "$SWAP_FILE"
    sudo chmod 600 "$SWAP_FILE"
    sudo mkswap "$SWAP_FILE" > /dev/null
    sudo swapon "$SWAP_FILE"
  else
    echo "::warning::Not enough free disk for a ${SWAP_GB} GB swap file; ISOCALL_ALIGN may be killed for lack of memory (exit 137)."
  fi
fi

if [[ "$ENGINE" == "docker" ]]; then
  GUARD_LOG="${RUNNER_TEMP}/docker_disk_guard.log"
  nohup bash "${SCRIPT_DIR}/ci_docker_disk_guard.sh" > "$GUARD_LOG" 2>&1 &
  echo "DOCKER_DISK_GUARD_LOG=${GUARD_LOG}" >> "$GITHUB_ENV"
  echo "Docker disk guard started (pid $!, log: ${GUARD_LOG}); Docker stores images in $(docker info --format '{{.DockerRootDir}}' 2> /dev/null || echo '?')"
fi

free -h
swapon --show
df -h / /mnt 2> /dev/null || df -h /
