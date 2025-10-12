python3 /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/prep_pf2_h5_splits.py \
  --csv /home/rover2/HW1_popgen/kursov/kursovaya/ML/Table1_for_ML.final_for_training.csv \
  --h5-root /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset \
  --outdir  /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass \
  --test-size 0.2


sudo rm -r /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result && pathogenfinder2 train -c /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/train.json -o /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result

python3 /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/extract_ytrue_ypred.py /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2/val_predictions.tsv


python3 /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/metrics_from_tsv.py \
  /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2/val_probs.tsv \
  --outdir /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2

```bash
sudo rm -r /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result && pathogenfinder2 train -c /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/train.json -o /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result \
&& python3 /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/extract_ytrue_ypred.py /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2/val_predictions.tsv \
&& python3 /home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/metrics_from_tsv.py \
  /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2/val_probs.tsv \
  --outdir /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2
  
```