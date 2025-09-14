# Step 1:

```bash
bash /home/rover2/HW1_popgen/kursov/PatogenPrediction/combine_fna.sh
```

# Step 2:

```bash
pathogenfinder2 predict \
  --multipleFiles "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna.txt" \
  -f genome \
  -o "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/result" \
  --prodigalPath /usr/local/bin/prodigal \
  --embedProteome map \
  --attProteins align \
  --dbProteins "/home/rover2/HW1_popgen/kursov/PatogenPrediction/uniref50.dmnd" \
  --diamondPath /home/rover2/HW1_popgen/kursov/PatogenPrediction/Diamond/diamond
```

# Step 2.5 (Optonoal to Test):
```bash
pathogenfinder2 predict \
  --multipleFiles "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna_test.txt" \
  -f genome \
  -o "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/result_test" \
  --prodigalPath /usr/local/bin/prodigal \
  --embedProteome map \
  --attProteins align \
  --dbProteins "/home/rover2/HW1_popgen/kursov/PatogenPrediction/uniref50.dmnd" \
  --diamondPath /home/rover2/HW1_popgen/kursov/PatogenPrediction/Diamond/diamond
```

Step 3:
```
visulalizations.ipynb
```