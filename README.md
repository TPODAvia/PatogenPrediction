# Step 1 (Optionally for Test):

```bash
pathogenfinder2 predict \
  -i /home/rover2/2856_genomes/prokka_outputs/195.363/195.363.faa \
  -f proteome \
  -o /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/result_test_single
```

or

```bash
pathogenfinder2 predict \
  --multipleFiles "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna_test.txt" \
  -f proteome \
  -o /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/result_test_multiple
```

# Step 2:

--dbProteins must point to a DIAMOND-indexed database (.dmnd), not a plain FASTA. Build it first, then pass the .dmnd path.

```bash
/home/rover2/HW1_popgen/kursov/PatogenPrediction/Diamond/diamond makedb \
  --in /home/rover2/HW1_popgen/kursov/PatogenPrediction/uniref50.fasta \
  -d /home/rover2/HW1_popgen/kursov/PatogenPrediction/uniref50
```

# Step 3:

```bash
bash /home/rover2/HW1_popgen/kursov/PatogenPrediction/combine_faa.sh
```

```bash
# Prevent Out of Memory
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True,max_split_size_mb:64,garbage_collection_threshold:0.8"
export TOKENIZERS_PARALLELISM=false

bash /home/rover2/HW1_popgen/kursov/PatogenPrediction/run_pf2_batches.sh
```

Step 4:
```
visulalizations.ipynb
```