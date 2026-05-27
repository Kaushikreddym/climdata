#!/usr/bin/env bash
# ============================================================================
# Download all ensemble members for GCFS22 / DWD-EPISODES2022 / svh20230401
# from the DWD ESGF THREDDS server.
#
# Variables: pr, hurs, tas, tasmax, tasmin, rsds   (frequency: day)
# Grid:      DE-0075x005   Cast package: svh20230401
#
# The script auto-detects:
#   - available sfc<DATE>  initialization sub-directories
#   - available r<N>i1p1   ensemble members
#   - available v<DATE>    version directories
#   - the exact filename(s) per (variable, ensemble, sfc, version)
# by parsing the THREDDS catalog.xml at each level.
#
# Requirements: bash, curl, grep, sed   (wget is used for the actual download
# because it resumes nicely; falls back to curl if wget is missing).
# ============================================================================

set -u
set -o pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HOST="https://esgf-data.dwd.de"
CATALOG_BASE="${HOST}/thredds/catalog/esgf3/data/climatepredictionsde/seasonal/output/public/DE-0075x005/DWD/GCFS22"
FILE_BASE="${HOST}/thredds/fileServer/esgf3/data/climatepredictionsde/seasonal/output/public/DE-0075x005/DWD/GCFS22"

CAST_PACKAGES=()
for year in {2023}; do
    for month in {01..12}; do
        CAST_PACKAGES+=("svh${year}${month}01")
    done
done
RCM="DWD-EPISODES2022"
RCM_VER="v1-r1"
FREQ="day"
#VARIABLES=("pr" "hurs" "tas" "tasmax" "tasmin" "rsds")
VARIABLES=("sfcWind")
OUTDIR="${OUTDIR:-/data01/FDS/muduchuru/Atmos/DWD/S2S/GCFS22_2010-2026}"
LOGFILE="${OUTDIR}/download.log"

# Optional: restrict to a subset of ensembles, e.g.  ENSEMBLES_FILTER="r1i1p1 r2i1p1"
ENSEMBLES_FILTER="${ENSEMBLES_FILTER:-}"
# Optional: dry run (list URLs only, do not download)
DRYRUN="${DRYRUN:-0}"
# Parallel downloads (1 = sequential, kind to the server)
PARALLEL="${PARALLEL:-1}"

mkdir -p "$OUTDIR"
: > "$LOGFILE"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGFILE" >&2; }
die()  { log "ERROR: $*"; exit 1; }

# Pick downloader
if   command -v wget >/dev/null 2>&1; then DL_TOOL="wget"
elif command -v curl >/dev/null 2>&1; then DL_TOOL="curl"
else die "Need either 'wget' or 'curl' on PATH."
fi
command -v curl >/dev/null 2>&1 || die "Need 'curl' for catalog scraping."

# Fetch a THREDDS catalog.xml and return child dataset names
# Args: <catalog-dir-url-without-trailing-slash>
list_children() {
    local url="$1/catalog.xml"
    # THREDDS catalogs use <dataset name="X" ID="..."> or <catalogRef xlink:title="X" .../>
    # We extract both 'name=' (datasets) and 'xlink:title=' (catalogRefs).
    curl -s -f "$url" 2>>"$LOGFILE" \
        | tr '\n' ' ' \
        | grep -oE '(xlink:title|name)="[^"]+"' \
        | sed -E 's/.*="([^"]+)"/\1/'
}

# List files (leaf datasets, .nc) in a catalog directory
list_files() {
    local url="$1/catalog.xml"
    curl -s -f "$url" 2>>"$LOGFILE" \
        | tr '\n' ' ' \
        | grep -oE 'urlPath="[^"]+\.nc"' \
        | sed -E 's/urlPath="([^"]+)"/\1/'
}

download_one() {
    local file_url="$1"
    local target="$2"
    if [[ "$DRYRUN" == "1" ]]; then
        echo "DRYRUN  $file_url"
        return 0
    fi
    if [[ -s "$target" ]]; then
        log "SKIP (exists) $(basename "$target")"
        return 0
    fi
    log "GET  $(basename "$target")"
    local tmp="${target}.part"
    if [[ "$DL_TOOL" == "wget" ]]; then
        wget -q --show-progress -c -O "$tmp" "$file_url" \
            && mv "$tmp" "$target" \
            || { rm -f "$tmp"; log "FAIL $file_url"; return 1; }
    else
        curl -L --fail --retry 3 --retry-delay 5 -C - -o "$tmp" "$file_url" \
            && mv "$tmp" "$target" \
            || { rm -f "$tmp"; log "FAIL $file_url"; return 1; }
    fi
}

# ---------------------------------------------------------------------------
# 1. Discover sfc<DATE> sub-directories under each CAST_PACKAGE
# ---------------------------------------------------------------------------
log "Scanning CAST_PACKAGES (${#CAST_PACKAGES[@]} packages for 2010-2026)..."

# ---------------------------------------------------------------------------
# 2. Loop over each CAST_PACKAGE -> sfc<DATE> -> ensembles -> versions per var
# ---------------------------------------------------------------------------
TOTAL_OK=0
TOTAL_FAIL=0

for CAST_PACKAGE in "${CAST_PACKAGES[@]}"; do
    log "==== CAST_PACKAGE: $CAST_PACKAGE ===="
    
    log "Scanning ${CAST_PACKAGE} for date directories matching *{2010-2026}{01-12}01..."
    SFC_DIRS=()
    # Discover all actual directories under CAST_PACKAGE
    while IFS= read -r entry; do
        # Match any prefix followed by YYYYMMDD01 pattern for years 2010-2026
        if [[ "$entry" =~ ^[a-z]{3}(201[0-9]|202[0-6])(0[1-9]|1[0-2])01$ ]]; then
            SFC_DIRS+=("$entry")
        fi
    done < <(list_children "${CATALOG_BASE}/${CAST_PACKAGE}" 2>/dev/null | sort -u)

    if [[ ${#SFC_DIRS[@]} -eq 0 ]]; then
        log "  No date directories found under ${CAST_PACKAGE} — skipping."
        continue
    fi
    log "Found ${#SFC_DIRS[@]} date dirs"

    # 2a. Loop over each sfc<DATE> -> discover ensembles -> discover versions per var
    for SFC in "${SFC_DIRS[@]}"; do
        log "---- sfc dir: $SFC ----"

        # 2b. Discover ensembles (rXiYpZ pattern)
        ENSEMBLES=()
        while IFS= read -r entry; do
            [[ "$entry" =~ ^r[0-9]+i[0-9]+p[0-9]+$ ]] && ENSEMBLES+=("$entry")
        done < <(list_children "${CATALOG_BASE}/${CAST_PACKAGE}/${SFC}" 2>/dev/null | sort -u)

        if [[ ${#ENSEMBLES[@]} -eq 0 ]]; then
            log "  No ensembles found for ${SFC} — skipping."
            continue
        fi

        # Optional filter
        if [[ -n "$ENSEMBLES_FILTER" ]]; then
            FILTERED=()
            for e in "${ENSEMBLES[@]}"; do
                for f in $ENSEMBLES_FILTER; do
                    [[ "$e" == "$f" ]] && FILTERED+=("$e")
                done
            done
            ENSEMBLES=("${FILTERED[@]}")
        fi

        log "  Ensembles (${#ENSEMBLES[@]}): ${ENSEMBLES[*]}"

        for ENS in "${ENSEMBLES[@]}"; do
            ENS_DIR_OUT="${OUTDIR}/${CAST_PACKAGE}/${SFC}/${ENS}"
            mkdir -p "$ENS_DIR_OUT"

            for VAR in "${VARIABLES[@]}"; do
                VAR_BASE="${CATALOG_BASE}/${CAST_PACKAGE}/${SFC}/${ENS}/${RCM}/${RCM_VER}/${FREQ}/${VAR}"

                # 2c. Discover version directories (v<DATE>) for this variable
                VERSIONS=()
                while IFS= read -r entry; do
                    [[ "$entry" =~ ^v[0-9]{8}$ ]] && VERSIONS+=("$entry")
                done < <(list_children "$VAR_BASE" 2>/dev/null | sort -u)

                if [[ ${#VERSIONS[@]} -eq 0 ]]; then
                    log "    [${ENS}/${VAR}] no version dir found — variable not published?"
                    continue
                fi

                # Use the latest version (lexicographic = chronological for vYYYYMMDD)
                LATEST_VER="${VERSIONS[${#VERSIONS[@]}-1]}"

                # 2d. List the actual .nc files under that version
                FILES=()
                while IFS= read -r f; do
                    FILES+=("$f")
                done < <(list_files "${VAR_BASE}/${LATEST_VER}")

                if [[ ${#FILES[@]} -eq 0 ]]; then
                    log "    [${ENS}/${VAR}/${LATEST_VER}] catalog empty — skipping."
                    continue
                fi

                for URLPATH in "${FILES[@]}"; do
                    FNAME="$(basename "$URLPATH")"
                    FILE_URL="${HOST}/thredds/fileServer/${URLPATH}"
                    TARGET="${ENS_DIR_OUT}/${FNAME}"
                    if download_one "$FILE_URL" "$TARGET"; then
                        TOTAL_OK=$((TOTAL_OK + 1))
                    else
                        TOTAL_FAIL=$((TOTAL_FAIL + 1))
                    fi
                done
            done
        done
    done
done

log "============================================================"
log "Done. Successful: ${TOTAL_OK}   Failed: ${TOTAL_FAIL}"
log "Output dir : ${OUTDIR}"
log "Log file   : ${LOGFILE}"
log "============================================================"

[[ "$TOTAL_FAIL" -eq 0 ]] || exit 2
exit 0
