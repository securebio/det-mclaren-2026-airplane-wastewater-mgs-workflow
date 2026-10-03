#!/usr/bin/env bash
# Build the SortMeRNA index that RUN's ribosomal screen (modules/local/sortmerna) uses,
# from SortMeRNA's standard sensitive database: SILVA SSU and LSU, and Rfam 5S and 5.8S
# (smr_v4.3_sensitive_db.fasta, from the v4.3.4 release's database.tar.gz).
#
# Usage:
#   bin/build_sortmerna_index.sh OUT_DIR [THREADS]
#
# Writes to OUT_DIR:
#   smr_v4.3_sensitive_db.fasta   the reference FASTA (params.sortmerna_ref)
#   idx/                          the index (params.sortmerna_index_dir)
#   index-info.txt                SortMeRNA version, container image, database release,
#                                 and checksums
#   build.log                     SortMeRNA's indexing log
# Upload OUT_DIR as a whole (e.g. aws s3 sync OUT_DIR s3://bucket/prefix/), then point
# params.sortmerna_ref and params.sortmerna_index_dir at the uploaded copies.
#
# The index is tied to the SortMeRNA version, so it is built with the container image the
# pipeline uses (the sortmerna label in configs/containers.config), via Docker. For the full
# database, building took about 30 minutes on 3-4 threads, and the index is about 4.0 GB.
#
# SortMeRNA names index files by a hash of the reference FASTA's file name, so the FASTA is
# mounted as ref.fasta for the build, and the pipeline links it under that name at run time.
#
# Options, as environment variables:
#   DB_TARBALL  a local copy of database.tar.gz, to skip the download (its checksum must
#               match the v4.3.4 release's)
#   FASTA       a FASTA to index instead of the release database, e.g. to regenerate the
#               test index from test-data/sortmerna/test_ref.fasta

set -euo pipefail

OUT_DIR="${1:?Usage: $0 OUT_DIR [THREADS]}"
THREADS="${2:-4}"
DB_RELEASE="v4.3.4"
DB_URL="https://github.com/sortmerna/sortmerna/releases/download/${DB_RELEASE}/database.tar.gz"
DB_SHA256="7bf684b9a0eb37e9b74c97e92f5a1f9acb9175f6fab03812d458da9049ce9660"
DB_FASTA_NAME="smr_v4.3_sensitive_db.fasta"
REF_PATH="/tmp/sortmerna-ref/ref.fasta"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# Print a message, and the tail of a log if given, and stop
fail() {
    echo "$1" >&2
    if [[ -f "${2:-}" ]]; then
        tail -n 20 "$2" >&2
    fi
    exit 1
}

# The image from the sortmerna block of configs/containers.config
IMAGE=$(awk '
    /withLabel:[[:space:]]*sortmerna[[:space:]]*\{/ { in_block = 1; next }
    in_block && /container[[:space:]]*=/ {
        if (match($0, /"[^"]+"/)) print substr($0, RSTART + 1, RLENGTH - 2)
        exit
    }
    in_block && /^[[:space:]]*\}/ { exit }
' "${REPO_ROOT}/configs/containers.config")
[[ -n "$IMAGE" ]] || fail "No sortmerna container in configs/containers.config"

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
if [[ -e "${OUT_DIR}/idx" ]]; then
    fail "${OUT_DIR}/idx already exists; remove it or choose another OUT_DIR"
fi
mkdir "${OUT_DIR}/idx"

# 1. Get the reference FASTA: a given one, or the sensitive database from the release
if [[ -n "${FASTA:-}" ]]; then
    fasta_name=$(basename "$FASTA")
    cp "$FASTA" "${OUT_DIR}/${fasta_name}"
    reference_source="${FASTA} (provided)"
    tarball_sha256="NA"
else
    fasta_name="$DB_FASTA_NAME"
    tarball="${DB_TARBALL:-${OUT_DIR}/database.tar.gz}"
    if [[ -z "${DB_TARBALL:-}" ]]; then
        echo "Downloading ${DB_URL}"
        curl -fsSL -o "$tarball" "$DB_URL"
    fi
    tarball_sha256=$(sha256sum "$tarball" | cut -d' ' -f1)
    if [[ "$tarball_sha256" != "$DB_SHA256" ]]; then
        fail "${tarball} has sha256 ${tarball_sha256}, not the ${DB_RELEASE} release's ${DB_SHA256}"
    fi
    fasta_member=$(tar -tzf "$tarball" | grep -E "(^|/)${DB_FASTA_NAME//./\\.}$")
    tar -xzf "$tarball" -O "$fasta_member" > "${OUT_DIR}/${fasta_name}"
    reference_source="SortMeRNA ${DB_RELEASE} database.tar.gz member ${fasta_member} (${DB_URL})"
    if [[ -z "${DB_TARBALL:-}" ]]; then
        rm "$tarball"
    fi
fi

# 2. Build the index in the pipeline's container
docker pull -q "$IMAGE" > /dev/null
echo "Building index with ${IMAGE} on ${THREADS} threads"
docker run --rm -u "$(id -u):$(id -g)" \
    -v "${OUT_DIR}/${fasta_name}:${REF_PATH}:ro" \
    -v "${OUT_DIR}/idx:/idx" \
    "$IMAGE" \
    sortmerna --ref "$REF_PATH" --idx-dir /idx --task 5 --workdir /tmp/sortmerna-work \
        --threads "$THREADS" > "${OUT_DIR}/build.log" 2>&1 \
    || fail "Index build failed:" "${OUT_DIR}/build.log"

# 3. Check that the index loads and finds a read cut from the first reference sequence
check_dir=$(mktemp -d)
awk '/^>/{if (seq) exit; next} {seq = seq $0} END{print "@check\n" substr(seq, 1, 150) "\n+\n" sprintf("%150s", "") }' \
    "${OUT_DIR}/${fasta_name}" | sed '4s/ /I/g' > "${check_dir}/check.fastq"
docker run --rm -u "$(id -u):$(id -g)" \
    -v "${OUT_DIR}/${fasta_name}:${REF_PATH}:ro" \
    -v "${OUT_DIR}/idx:/idx:ro" \
    -v "${check_dir}:/check" \
    "$IMAGE" \
    sortmerna --ref "$REF_PATH" --idx-dir /idx --index 0 --reads /check/check.fastq \
        --fastx --aligned /check/aligned --other /check/other \
        --workdir /tmp/sortmerna-work --threads 1 > "${check_dir}/check.log" 2>&1 \
    || fail "Index check failed:" "${check_dir}/check.log"
# SortMeRNA names FASTQ output .fq
n_aligned=$(( $(wc -l < "${check_dir}/aligned.fq") / 4 ))
rm -rf "$check_dir"
[[ "$n_aligned" -eq 1 ]] || fail "Index check failed: the test read did not align"

# 4. Record what was built
version=$(docker run --rm "$IMAGE" sortmerna --version 2>&1 | sed -n 's/^SortMeRNA version //p')
digest=$(docker image inspect "$IMAGE" --format '{{index .RepoDigests 0}}')
cat > "${OUT_DIR}/index-info.txt" <<EOF_INFO
sortmerna_version: ${version}
container_image: ${IMAGE}
container_digest: ${digest}
reference_source: ${reference_source}
database_tarball_sha256: ${tarball_sha256}
reference_fasta: ${fasta_name}
reference_fasta_sha256: $(sha256sum "${OUT_DIR}/${fasta_name}" | cut -d' ' -f1)
reference_path_at_build: ${REF_PATH}
index_files: $(ls "${OUT_DIR}/idx" | wc -l) files, $(du -sb "${OUT_DIR}/idx" | cut -f1) bytes
built_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)
build_threads: ${THREADS}
EOF_INFO
cat "${OUT_DIR}/index-info.txt"
