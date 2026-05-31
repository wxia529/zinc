#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/fetch_deps.sh [--force] [OUT_DIR]

Download and unpack zinc's CMake FetchContent dependencies for offline builds.

Run this on a machine with network access, then copy OUT_DIR to the offline
machine together with the zinc source checkout.

Examples:
  scripts/fetch_deps.sh
  scripts/fetch_deps.sh /tmp/zinc-deps

Offline configure:
  cmake -C /path/to/zinc-deps/zinc-deps.cmake -S . -B build
USAGE
}

force=0
out_dir=".zinc-deps"

while (($#)); do
  case "$1" in
    --force)
      force=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      out_dir="$1"
      shift
      ;;
  esac
done

out_parent="$(dirname "$out_dir")"
mkdir -p "$out_parent"
out_dir="$(cd "$out_parent" && pwd)/$(basename "$out_dir")"
archive_dir="$out_dir/archives"
source_dir="$out_dir/src"
mkdir -p "$archive_dir" "$source_dir"

deps=(
  "eigen3|eigen-3.4.0.tar.gz|https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz|8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72|FETCHCONTENT_SOURCE_DIR_EIGEN3"
  "pugixml|pugixml-1.14.tar.gz|https://github.com/zeux/pugixml/archive/refs/tags/v1.14.tar.gz|610f98375424b5614754a6f34a491adbddaaec074e9044577d965160ec103d2e|FETCHCONTENT_SOURCE_DIR_PUGIXML"
  "CLI11|CLI11-2.4.2.tar.gz|https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.4.2.tar.gz|f2d893a65c3b1324c50d4e682c0cdc021dd0477ae2c048544f39eed6654b699a|FETCHCONTENT_SOURCE_DIR_CLI11"
  "fmt|fmt-10.2.1.tar.gz|https://github.com/fmtlib/fmt/archive/refs/tags/10.2.1.tar.gz|1250e4cc58bf06ee631567523f48848dc4596133e163f02615c97f78bab6c811|FETCHCONTENT_SOURCE_DIR_FMT"
  "googletest|googletest-1.14.0.tar.gz|https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz|8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7|FETCHCONTENT_SOURCE_DIR_GOOGLETEST"
)

download() {
  local url="$1"
  local path="$2"

  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --output "$path" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget --tries=3 --output-document="$path" "$url"
  else
    echo "Neither curl nor wget is available." >&2
    exit 1
  fi
}

verify_sha256() {
  local path="$1"
  local expected="$2"
  local actual

  if command -v sha256sum >/dev/null 2>&1; then
    actual="$(sha256sum "$path" | awk '{print $1}')"
  elif command -v shasum >/dev/null 2>&1; then
    actual="$(shasum -a 256 "$path" | awk '{print $1}')"
  else
    echo "Neither sha256sum nor shasum is available." >&2
    exit 1
  fi

  if [[ "$actual" != "$expected" ]]; then
    echo "SHA256 mismatch for $path" >&2
    echo "  expected: $expected" >&2
    echo "  actual:   $actual" >&2
    exit 1
  fi
}

extract_archive() {
  local archive="$1"
  local target="$2"
  local tmp_dir
  local top_dir

  tmp_dir="$(mktemp -d)"
  cmake -E chdir "$tmp_dir" cmake -E tar xzf "$archive"
  top_dir="$(find "$tmp_dir" -mindepth 1 -maxdepth 1 -type d | head -n 1)"

  if [[ -z "$top_dir" ]]; then
    echo "Could not find extracted top-level directory in $archive" >&2
    rm -rf "$tmp_dir"
    exit 1
  fi

  rm -rf "$target"
  mkdir -p "$(dirname "$target")"
  mv "$top_dir" "$target"
  rm -rf "$tmp_dir"
}

for dep in "${deps[@]}"; do
  IFS='|' read -r name archive_name url sha256 cmake_var <<<"$dep"
  archive="$archive_dir/$archive_name"
  target="$source_dir/$name"

  if [[ "$force" -eq 1 ]]; then
    rm -f "$archive"
    rm -rf "$target"
  fi

  if [[ ! -f "$archive" ]]; then
    echo "Downloading $name"
    download "$url" "$archive"
  else
    echo "Using existing archive for $name"
  fi

  verify_sha256 "$archive" "$sha256"

  if [[ ! -d "$target" ]]; then
    echo "Extracting $name"
    extract_archive "$archive" "$target"
  else
    echo "Using existing source directory for $name"
  fi
done

cache_file="$out_dir/zinc-deps.cmake"
cat >"$cache_file" <<'CMAKE'
get_filename_component(_ZINC_DEPS_DIR "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

set(FETCHCONTENT_FULLY_DISCONNECTED ON CACHE BOOL "Use pre-fetched zinc dependencies only" FORCE)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "Disable FetchContent dependency updates" FORCE)
set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON CACHE BOOL "Do not write user package registry entries" FORCE)

set(FETCHCONTENT_SOURCE_DIR_EIGEN3 "${_ZINC_DEPS_DIR}/src/eigen3" CACHE PATH "Offline Eigen source" FORCE)
set(FETCHCONTENT_SOURCE_DIR_PUGIXML "${_ZINC_DEPS_DIR}/src/pugixml" CACHE PATH "Offline pugixml source" FORCE)
set(FETCHCONTENT_SOURCE_DIR_CLI11 "${_ZINC_DEPS_DIR}/src/CLI11" CACHE PATH "Offline CLI11 source" FORCE)
set(FETCHCONTENT_SOURCE_DIR_FMT "${_ZINC_DEPS_DIR}/src/fmt" CACHE PATH "Offline fmt source" FORCE)
set(FETCHCONTENT_SOURCE_DIR_GOOGLETEST "${_ZINC_DEPS_DIR}/src/googletest" CACHE PATH "Offline googletest source" FORCE)

message(STATUS "Using zinc offline dependencies from ${_ZINC_DEPS_DIR}")
CMAKE

cat <<EOF

Offline dependency bundle is ready:
  $out_dir

Copy that directory to the offline machine and configure zinc with:
  cmake -C $out_dir/zinc-deps.cmake -S . -B build

EOF
