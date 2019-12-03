#Dependencies:
- Python 3.7
- matplotlib
- numpy
- pandas
- https://github.com/breecummins/min_interval_posets.git


#Analyses
This code is to repeat the partial ordering analysis of the malaria strain and mouse tissue data in An Intrinsic Oscillator Drives the Blood Stage Cycle of the Malaria Parasite, Plasmodium falciparum, Smith et al. (submitted Dec 2019).

This computation is cpu, memory, and time intensive. It's recommended to be performed on a cluster.


##To create the data and supplemental figures for the malaria strains, run

```python malaria_scripts.py define_params_mal <number of cores>```

A date and time-stamped results folder is the output. 

To make the values in the supplemental table, change directories into the results folder and run

`python ../analyze.py`



##To create the data and supplemental figures for the mouse tissues, run

`python malaria_scripts.py define_params_mouse <number of cores>`

A date and time-stamped results folder is the output. 

To make the values in the supplemental table, go to `analyze.py` and comment out the line
`load_files(strains = ["D6","FVO","SA250"])`
and uncomment the line
`# load_files(strains = ["lung","kidney"])`
Then change directories into the results folder and run

`python ../analyze.py`





#Data files
Analyses performed by Kim Roche and Lauren Smith

####WRAPPED DOWNSAMPLED EXPRESSION TIME SERIES DATA FOR P. falciparum:

downsampled_{strain}_expression_wrapped_offset.csv
                         Time-series RNA-seq in FPKM. Gene expression interpolated via PCHIP and wrapped to a lowest-error single cycle. Offset to a common starting point (~%50/50 troph-schizont transition)

####TRUNCATED EXPRESSION TIME SERIES FOR MOUSE TISSUE DATA:

{tissue}_in_phase_0_05.tsv


####GENELIST FOR MALARIA:

allstrains_p25_genelist.txt
                         Set of genes periodic in all strains. The "baseline" set of genes.

From this baseline set, peak time was transformed to "percent of cycle completed" based on wrapped expression data/wrapped staging data. 

3D7_allstrain_outofphase_0.05_two.txt
			Genes differing by < 5% of the cycle in at least TWO strains compared to 3D7. Filters for minimum mis-ordering reduce gene set size.


####GENELIST FOR MOUSE:

mouse_gene_in_phase.txt
            Genes differing by < 5% of the cycle in at least one tissue compared to liver. Filters for minimum mis-ordering reduce gene set size.
            
            
#Shifted data files

All files ending in `_shifted.tsv` were produced by the script `make_phase_shift.py` to define the baseline comparison in the manuscript. Rerunning

`python make_phase_shift.py <data file name>`

will create a new baseline for that data file.