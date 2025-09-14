# write absolute paths into this list
LIST="/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna.txt"
mkdir -p /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset

find /home/rover2/2856_genomes/mcSEED_2856_fna \
  -maxdepth 1 -type f \( -iname '*.fna' -o -iname '*.fa' -o -iname '*.fasta' \) \
  -print0 | sort -z | tr '\0' '\n' > "$LIST"

wc -l "$LIST"; head -n 3 "$LIST"
