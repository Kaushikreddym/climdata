#!/bin/bash

start_year=""
end_year=""
root_dir=""
parallel_jobs=4
data_type="daily"  # default

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --year_start) start_year="$2"; shift 2 ;;
        --year_end) end_year="$2"; shift 2 ;;
        --path) root_dir="$2"; shift 2 ;;
        --parallel) parallel_jobs="$2"; shift 2 ;;
        --type) data_type="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate
if [[ -z "$start_year" || -z "$end_year" || -z "$root_dir" ]]; then
    echo "Usage: $0 --year_start <year> --year_end <year> --path <root_dir> [--parallel N] [--type daily|hourly]"
    exit 1
fi

# Base URLs
if [[ "$data_type" == "daily" ]]; then
    base_url="https://www.ncei.noaa.gov/data/global-summary-of-the-day/archive"
elif [[ "$data_type" == "hourly" ]]; then
    base_url="https://www.ncei.noaa.gov/data/global-hourly/archive/csv"
else
    echo "Invalid --type, must be 'daily' or 'hourly'"
    exit 1
fi

archive_dir="${root_dir}/archive"
mkdir -p "$archive_dir"

download_and_extract() {
    year="$1"
    file="${archive_dir}/${year}.tar.gz"
    extract_dir="${root_dir}/${year}"

    url="${base_url}/${year}.tar.gz"
    echo "--> Processing $data_type year: $year"

    # Download only if missing
    if [[ ! -f "$file" ]]; then
        echo "Downloading: $url"
        wget -c -O "$file" "$url"
    else
        echo "Archive exists: $file (skip download)"
    fi

    # Create extraction folder
    mkdir -p "$extract_dir"

    # Skip extraction if folder has contents
    if [[ "$(ls -A "$extract_dir")" ]]; then
        echo "Already extracted: $extract_dir (skip)"
    else
        echo "Extracting CSV into: $extract_dir ..."
        tar -xzf "$file" -C "$extract_dir" --wildcards "*.csv"
    fi
}

export -f download_and_extract
export base_url archive_dir root_dir data_type

# Parallel execution if possible
if command -v parallel > /dev/null 2>&1; then
    seq "$start_year" "$end_year" | parallel -j "$parallel_jobs" download_and_extract {}
else
    jobs_count=0
    for year in $(seq "$start_year" "$end_year"); do
        download_and_extract "$year" &
        ((jobs_count++))
        if (( jobs_count >= parallel_jobs )); then
            wait -n
            ((jobs_count--))
        fi
    done
    wait
fi

echo "✅ Completed ✅"
echo "📦 Archives: $archive_dir"
echo "📂 Extracted Data: $root_dir/<year>"

