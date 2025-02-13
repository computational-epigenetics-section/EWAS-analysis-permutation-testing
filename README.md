# Permutation Testing Suite

This is a suite which enables to calculate the enrichment ratios of 100,000 permuted sets of CoRSIVs and control probes.

## File Structure

- `Data` - This has all the input files for the permutation testing including the probe lists for each disease.

- `permuted_sets` - This includes 100,000 files. Each files contains permuted CoRSIV and control probes. This is the output of `get_permutations.py`. Due to the large size of this file, it was not added here. It's located in the SPHERE cluster (/home/maduranga/Permutation_testing_code_SPHERE).

- `p1.py` to `p10.py` - Calculate enrichment ratios for 100,000 permuted sets. Each file processes 10,000 sets which can be run in parallel if needed.

- `merge_chunks.py` - Concatenates the chunks to create one final output for each disease.



## Usage

1) Add the probe list CSV file of the disease to `Data` folder.

- The file name should end with "_all_probes.csv". Ex: `neurological_all_probes.csv` or `Insulin Resistance_all_probes.csv` or `Diabetes Mellitus, Type 2_all_probes.csv`

2) Use the `permutation_testing.slrum` to run permutation testing.

- Example:

```
#!/bin/sh
#SBATCH --time=48:00:00 -n24 -p bynode

for d in "neurological"
do
    python p1.py "${d}"
    python p2.py "${d}"
    python p3.py "${d}"
    python p4.py "${d}"
    python p5.py "${d}"
    python p6.py "${d}"
    python p7.py "${d}"
    python p8.py "${d}"
    python p9.py "${d}"
    python p10.py "${d}"
done
```

### OR

You can run each of `p1.py` to `p10.py` in parallel.

- Example:

```
#!/bin/sh
#SBATCH --time=24:00:00 -n24 -p bynode

python p1.py "neurological"
   
```

```
#!/bin/sh
#SBATCH --time=24:00:00 -n24 -p bynode

python p2.py "neurological"
   
```
###   .  .  .  .  .  .  .  .  .  .  .  .  
###   .  .  .  .  .  .  .  .  .  .  .  .  

###   .  .  .  .  .  .  .  .  .  .  .  .  




```
#!/bin/sh
#SBATCH --time=24:00:00 -n24 -p bynode

python p10.py "neurological"
   
```

3) After the above script finishes running, merge the results by running `merge_chunks.slrum`.

- Example:

```
#!/bin/sh
#SBATCH --time=48:00:00 -n24 -p bynode

python merge_chuncks.py "neurological"

```

The output will be generated in `final_permutation_testing_results` folder.
