#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: wget_frags.sh <job_id>"
  exit 1
fi

# This script download fragments from the robetta server.  Really, this script 
# just exists so that you don't have to remember which options you need to pass 
# to wget to have it to download a whole directory.

wget \
  --recursive \
  --no-parent \
  --reject 'index*' \
  --execute robots=off \
  http://old.robetta.org/downloads/fragments/$1/
