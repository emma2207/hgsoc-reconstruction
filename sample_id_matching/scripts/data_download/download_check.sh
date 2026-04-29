#!/bin/bash

set -u

usage() {
	echo "Usage: $0 <files_directory> [checksum_file.txt]"
	echo "  files_directory    Directory containing files to verify"
	echo "  checksum_file.txt  Optional checksum list file (default: first md5_*.txt in current directory)"
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
	usage
	exit 2
fi

FILES_DIR="$1"
CHECKSUM_FILE="${2:-}"

if [[ ! -d "$FILES_DIR" ]]; then
	echo "ERROR: Directory not found: $FILES_DIR"
	exit 2
fi

if ! command -v md5sum >/dev/null 2>&1; then
	echo "ERROR: md5sum command not found in PATH."
	echo "Install coreutils or ensure md5sum is available."
	exit 2
fi

if [[ -z "$CHECKSUM_FILE" ]]; then
	# Use first checksum file matching md5_*.txt in current working directory.
	CHECKSUM_FILE="$(find . -maxdepth 1 -type f -name 'md5_*.txt' | head -n 1)"
fi

if [[ -z "$CHECKSUM_FILE" || ! -f "$CHECKSUM_FILE" ]]; then
	echo "ERROR: Could not find checksum file."
	echo "Pass it explicitly as the second argument."
	usage
	exit 2
fi

echo "Verifying checksums"
echo "  Files directory: $FILES_DIR"
echo "  Checksum file : $CHECKSUM_FILE"
echo

total=0
passed=0
failed=0
missing=0

while read -r expected_md5 file_name; do
	[[ -z "${expected_md5:-}" ]] && continue
	[[ "${expected_md5:0:1}" == "#" ]] && continue

	total=$((total + 1))
	file_path="$FILES_DIR/$file_name"

	if [[ ! -f "$file_path" ]]; then
		echo "WARNING: Missing file -> $file_name"
		missing=$((missing + 1))
		failed=$((failed + 1))
		continue
	fi

	actual_md5="$(md5sum "$file_path" | awk '{print $1}')"

	if [[ "$actual_md5" == "$expected_md5" ]]; then
		echo "OK: $file_name"
		passed=$((passed + 1))
	else
		echo "WARNING: Checksum mismatch -> $file_name"
		echo "  expected: $expected_md5"
		echo "  actual  : $actual_md5"
		failed=$((failed + 1))
	fi
done < "$CHECKSUM_FILE"

echo
echo "Summary:"
echo "  Total checked : $total"
echo "  Passed        : $passed"
echo "  Failed        : $failed"
echo "  Missing files : $missing"

if [[ "$failed" -gt 0 ]]; then
	echo "WARNING: One or more files failed checksum verification."
	exit 1
fi

echo "All files passed checksum verification."
exit 0

