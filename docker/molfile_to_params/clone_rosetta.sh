#!/usr/bin/env bash
set -euo pipefail

target_dir=${1:-}
tag=${2:-}

if [ -z "$target_dir" -o -z "$tag" ]; then
  echo "Usage: ./clone_rosetta.sh <target_dir> <tag>"
  exit 1
fi

# https://stackoverflow.com/questions/600079/how-do-i-clone-a-subdirectory-only-of-a-git-repository
git clone \
  --depth 1 \
  --branch $tag \
  --config advice.detachedHead=0 \
  --filter=blob:none \
  --sparse \
  git@github.com:RosettaCommons/main.git \
  $target_dir

cd $target_dir
git sparse-checkout set 'source/scripts/python/public'

 
