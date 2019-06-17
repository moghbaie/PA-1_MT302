# Line 1 PA-1 vs. MT302

##### contact moghbaie@rockefeller.edu for any questions
##### 06/17/2019


Project code and data is organized in following order:

* [src](https://github.com/moghbaie/L1_CRC_v2/tree/master/src)
* [Input_data](https://github.com/moghbaie/L1_CRC_v2/tree/master/Input_data/ReadMe_Input.md)
* [Image](https://github.com/moghbaie/L1_CRC_v2/tree/master/Image)


## Project Pipeline:
<img src="https://github.com/moghbaie/L1_CRC_IP_MS/blob/master/CRC_pipeline.png" alt="Pipeline" width="1000"></img>


## QC:

## Imputation:
### Comparison before - after imputation:
<img src="Image/Imputation/comparison_before_after_imputation_Colon_6_2.png" alt="Pipeline" width="900"></img>
<br>
<img src="Image/Imputation/comparison_before_after_imputation_Liver_6_2.png" alt="Pipeline" width="900"></img>
<br>
<img src="Image/Imputation/comparison_before_after_imputation_Ovary_6_2.png" alt="Pipeline" width="600"></img>

### Average intensities before - after imputation:
![Average intensities before imputation](Image/Imputation/average_intensities_before_imputation.png)
![Average intensities after imputation](Image/Imputation/average_intensities_after_imputation.png)

## ANOVA result: 

### Venndiagram
![VennDiagram of significant proteins ](Image/Volcano_plot/Significant_VennDiagram.png)

### Volcano plots and mutations
#### comparison between Tumor C (colon) against ORF1 and IgG
![Volcano plot compare Colon Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_163T_IgG_7.png)
<br>
#### comparison between Tumor C (colon) and normal tissue against ORF1 
![Volcano plot compare Colon Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_163N_ORF1_8.png)
<br>
#### comparison between Tumor B (liver) against ORF1 and IgG 
![Volcano plot compare Liver Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_159T_IgG_4.png)
##### Mutated peptides in significant ORF1 loci in tumor B against ORF1 and IgG 
![Observed significant Liver colon tissue against ORF1 and IgG ](Image/Phylogenetic_tree/Liver_Tumor_159T_IgG_4.png)
##### Mutation frequency in tumor B  against ORF1 and IgG
![Mutation frequency in liver tumor  against ORF1 and IgG ](Image/Phylogenetic_tree/mutation_fre_Liver_Tumor_159T_IgG_4.png)
<br>
#### comparison between Tumor B (liver) and normal tissue against ORF1
![Volcano plot compare Liver Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_159N_ORF1_5.png)
##### Mutated peptides in significant ORF1 loci in Tumor B and normal tissue against ORF1
![Observed significant Liver Tumor tissue and normal tissue against ORF1 ](Image/Phylogenetic_tree/Liver_Tumor_159N_ORF1_5.png)
##### Mutation frequency in tumor B and normal tissue against ORF1
![Mutation frequency in liver tumor and normal tissue against ORF1 ](Image/Phylogenetic_tree/mutation_fre_Liver_Tumor_159N_ORF1_5.png)
<br>
#### comparison between Tumor A (ovary) against ORF1 and IgG - phase 1
![Volcano plot compare Ovary Tumor tissue to Igg - phase1](Image/Volcano_plot/Volcano_plot_Ovary_Tumor_144T_IgG_2.png)
##### Mutated peptides in significant ORF1 loci in Tumor A against ORF1 and IgG - phase1
![Observed significant Ovary Tumor tissue against ORF1 and IgG - phase1](Image/Phylogenetic_tree/Ovary_Tumor_144T_IgG_2.png)
##### Mutation frequency in tumor A against ORF1 and IgG
![Mutation frequency in Ovary tumor  against ORF1 and IgG ](Image/Phylogenetic_tree/mutation_fre_Ovary_Tumor_144T_IgG_2.png)
<br>
#### comparison between Tumor A (ovary) against ORF1 and IgG - phase 2
![Volcano plot compare Ovary Tumor tissue to Igg - phase2](Image/Volcano_plot/Volcano_plot_Ovary_Tumor_144T_IgG_10.png)
##### Mutated peptides in significant ORF1 loci in Tumor A tissue against ORF1 and IgG - phase2
![Observed significant Ovary Tumor tissue against ORF1 and IgG - phase2](Image/Phylogenetic_tree/Ovary_Tumor_144T_IgG_10.png)
##### Mutation frequency in tumor A and normal tissue against ORF1
![Mutation frequency in Ovary tumor  and normal tissue against ORF1  ](Image/Phylogenetic_tree/mutation_fre_Ovary_Tumor_144T_IgG_10.png)

## Integrated plot:
### Venndiagram:
![VenDiagram of expressed proteins ](Image/Integrated_plot/VennDiagram.png)

### MDS plot:
![MDS plot ](Image/Integrated_plot/movie.gif)
<br>
<br>
### Polar map of proteins that were significant in two conditions or were observed in previous researches and were significant in at least one condition:
![Polar map - comparison of Colon significant proteins in different tissue](Image/Integrated_plot/heatmap_significant_eLife.png)

### Heat map of all significant proteins across different tissues:
![Heatmap - all significant proteins](Image/Integrated_plot/heatmap_all_significant.png)
"# PA-1_MT302" 
