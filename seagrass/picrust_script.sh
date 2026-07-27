#!/bin/bash

conda activate picrust2

less sequences.fna
biom head -i biom_table.biom
biom summarize-table -i biom_table.biom

mkdir output
cd output

# Place reads into reference tree 
# Change -p 6 to the number of cores to use

place_seqs.py -s ../sequences.fna -o out.tre -p 6 \
              --intermediate intermediate/place_seqs

<< 'COMMENT'
Warning - 42 input sequences aligned poorly to reference sequences (--min_align option specified a minimum proportion of 0.8 aligning to reference sequences). These input sequences will not be placed and will be excluded from downstream steps.
This is the set of poorly aligned input sequences to be excluded: ASV4125, ASV8109, ASV5629, ASV6929, ASV4998, ASV3676, ASV3223, ASV6871, ASV6782, ASV1552, ASV7809, ASV5055, ASV8025, ASV8053, ASV6474, ASV5989, ASV6950, ASV3705, ASV3523, ASV6593, ASV3624, ASV6478, ASV7566, ASV5181, ASV5862, ASV3882, ASV7629, ASV5172, ASV7946, ASV6036, ASV7153, ASV8026, ASV1149, ASV5234, ASV3696, ASV1697, ASV8174, ASV7947, ASV5876, ASV3909, ASV4090, ASV8047
COMMENT

# Hidden state prediction of gene families
# Change -p 6 to the number of cores to use
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 6 -n
hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 6
hsp.py -i KO -t out.tre -o KO_predicted.tsv.gz -p 6

<< 'COMMENT'
The script hsp.py can also output the nearest-sequenced taxon index (NSTI) values for each ASV (specified by the -n option), which correspond to the branch length in the tree from the placed ASV to the nearest reference 16S sequence. This metric is a rough guide for how similar an ASV is to an existing reference sequence. Predictions are more accurate for ASVs that have well-characterized close relatives.
ASVs with a NSTI score above 2 are usually noise. It can be useful to take a look at the distribution of NSTI values for your ASVs to determine how well-characterized your community is overall and whether there are any outliers.
COMMENT

zless -S marker_predicted_and_nsti.tsv.gz
zless -S EC_predicted.tsv.gz
zless -S KO_predicted.tsv.gz

# Generate metagenome predictions
<< 'COMMENT'
The read depth per ASV is divided by the predicted 16S copy numbers. This is performed to help control for variation in 16S copy numbers across organisms, which can result in interpretation issues. For instance, imagine an organism with five identical copies of the 16S gene that is at the same absolute abundance as an organism with one 16S gene. The ASV corresponding to the first organism would erroneously be inferred to be at higher relative abundance simply because this organism had more copies of the 16S gene. The ASV read depths per sample (after normalizing by 16S copy number) are multiplied by the predicted gene family copy numbers per ASV.
COMMENT

metagenome_pipeline.py -i ../biom_table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out --strat_out
metagenome_pipeline.py -i ../biom_table.biom -m marker_predicted_and_nsti.tsv.gz -f KO_predicted.tsv.gz \
                       -o KO_metagenome_out --strat_out
# 176 of 8155 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.

<< 'COMMENT'
EC_metagenome_out/weighted_nsti.tsv.gz - the mean NSTI value per sample (when taking into account the relative abundance of the ASVs). This file can be useful for identifying outlier samples in your dataset. In PICRUSt1 weighted NSTI values < 0.06 and > 0.15 were suggested as good and high, respectively. The cut-offs can be useful for getting a ball-park of how your samples compare to other datasets, but a weighted NSTI score > 0.15 does not necessarily mean that the predictions are meaningless.
COMMENT

zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
zless -S EC_metagenome_out/pred_metagenome_contrib.tsv.gz
zless -S KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
zless -S KO_metagenome_out/pred_metagenome_contrib.tsv.gz

<< 'COMMENT'
One thing to watch out for is that if you want to input this table of ASV contributions to another tool you should not use the --min_reads and --min_samples option to collapse rare ASVs to the RARE category since treating rare taxa as a single taxon would likely be inappropriate for many analyses (this will not happen by default).
COMMENT

# Pathway level inference
# Change -p 6 to the number of cores available for use
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o pathways_out -p 6
                    
# Functional descriptions 

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
