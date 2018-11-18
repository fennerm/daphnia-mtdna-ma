#!/usr/bin/env bash
# Extract sequences for single copy genes from the nuclear reference sequence.
set -euo pipefail

# TATA binding protein, X-box binding protein, CAPON/NOS1AP
TARGETS="APZ42_026392 APZ42_015309 APZ42_021922"

get_coordinates_from_gff() {
    local _gene_id="$1"
    local _gff_file="$2"
    # Second awk command removes duplicate lines
    grep "$_gene_id" "$_gff_file" | awk '{print $1, $4, $5}' | awk '!seen[$0]++'
}

extract_fasta_positions() {
    local _reference="$1"
    local _seq_id="$2"
    local _start="$3"
    local _stop="$4"
    echo "$(grep "$_seq_id" "$_reference")"
    extract_seq_from_fasta "$_seq_id" "$_reference" \
        | head -n -1 \
        | unwrap_fasta \
        | tail -n -1 \
        | cut -c "$_start"-"$_stop"
}

main() {
    gff="$1"
    reference="$2"
    
    for target in $TARGETS; do
        position="$(get_coordinates_from_gff "$target" "$gff")"
        seq_id="$(echo "$position" | awk '{print $1}')"
        start_pos="$(echo "$position" | awk '{print $2}')"
        stop_pos="$(echo "$position" | awk '{print $3}')"
        extract_fasta_positions "$reference" "$seq_id" "$start_pos" "$stop_pos"
    done
}


main "$@"
