<!--  Provide a title.         -->
## Major depression and atherosclerotic disease: Linking shared genetics to pathways in blood, brain, heart, and atherosclerotic plaques

<!--  Provide details on the people involved and the project ID.         -->
*Collaborators*

Emma Pruin, Meike Bartels, Ernest Diez Benavente, Noortje van den Dungen, Joost Hoekstra, Dominique D.P. de Kleijn, Michal Mokry, Gerard Pasterkamp, Brenda Penninx, Wouter Peyrot, Hester M. den Ruijter, S.W. van der Laan, PhD*, Y. Milaneschi, PhD*.
* these authors contributed equally

**Project ID** [`MR_CVD_MDD`]


### Background
<!--  Provide some background, study design, results, etc.         -->
The increased risk of atherosclerotic diseases (stroke, coronary artery disease [CAD]) observed in depression may stem from shared pathophysiology. We examined whether: 1) major depression (MD) and atherosclerotic traits share genetic risk, and 2) altered gene expression in various tissues linked to shared genetics has a potential causal role in depression etiology.

### Methods
Data from the largest genome-wide association studies of MD (N=3,887,532) and 8 atherosclerotic traits (N=26,909-1,308,460) were used in Mendelian randomization and colocalization to detect cross-trait causal associations and genomic loci containing shared causal variants. In shared loci, summary data-based Mendelian randomization estimated the effects of gene expression on MD etiology using expression quantitative trait loci datasets from whole blood, brain and heart tissues and atherosclerotic plaques from the Athero-Express Biobank Study.

### Results
MD genetic liability increased risk of any stroke (OR=1.15, p=9.47×10−8), ischemic stroke (OR=1.16, p=1.52×10−7), small vessel disease (OR=1.34, p=4.76×10−5) and CAD (OR=1.2, 95%CIs=1.13-1.26, p=3.76×10−22). Eight genomic regions harbored potentially shared causal variants, including one on chromosome 7 linking MD with any stroke, ischemic stroke and CAD. Altered expression of 16 genes in blood, 10 in brain, and 6 in heart was found causal for MD etiology. In atherosclerotic plaques, one gene was linked to MD at nominal significance only.

### Conclusion
Major depression and atherosclerotic diseases share genetic risk potentially acting in depression pathophysiology through expression of genes in blood, brain and heart tissues. An involvement of atherosclerotic plaques in depression etiology was not supported. Identified pathways could guide the development of new treatments to prevent depression-heightened atherosclerotic risk.


**Figure 1: Genotyped individuals in the Athero-Express Biobank Study**
![Genotyped individuals in the Athero-Express Biobank Study](PLOTS/20250328.overlap.AEDB_AEGS123.UpSetR.png)


#### Whole-plaque RNAseq

For the expression analysis we used carotid plaque-derived bulk RNAseq data and queried it for the gene list. Below a graph showing the overall expression of the genes (not all are in the data) compared to the mean expression of 1,000 randomly picked genes. 

**Figure 2: Overall expression of target genes in carotid plaques from the Athero-Express Biobank Study**
![Overall expression of target genes in carotid plaques from the Athero-Express Biobank Study](PLOTS/20250328.TargetExpression_vs_1000genes.png)


#### Single cell RNAseq

We projected target genes to the single-cell RNAseq data derived from 37 carotid plaque samples. We identified cell communities (Figure 2), mapped and projected target gene expression to the cell communities (Figure 3). 

**Figure 3: Cell communities identified in carotid plaques from the Athero-Express Biobank Study**
![Cell communities identified in carotid plaques from the Athero-Express Biobank Study](PLOTS/20250328.UMAP.png)


**Figure 4: Dotplot showing expression of target genes per cell type in carotid plaques from the Athero-Express Biobank Study**
![Dotplot showing expression of target genes per cell type in carotid plaques from the Athero-Express Biobank Study](PLOTS/20250328.DotPlot.Targets.png)


### Where do I start?

You can load this project in RStudio by opening the file called 'MR_CVD_MDD.Rproj'.

### Project structure

<!--  You can add rows to this table, using "|" to separate columns.         -->

File                                    | Description                          | Usage         
--------------------------------------- | ------------------------------------ | --------------
README.md                               | Description of project               | Human editable
MR_CVD_MDD.Rproj                        | Project file                         | Loads project
LICENSE                                 | User permissions                     | Read only
.worcs                                  | WORCS metadata YAML                  | Read only
renv.lock                               | Reproducible R environment           | Read only
images                                  | image directory for project          | Human editable
BASELINE                                | Baseline characteristics directory   | Human editable
OUTPUT                                  | Output directory                     | Human editable
PLOTS                                   | Some results                         | Human editable
scripts                                 | Scripts directory                    | Human editable
targets                                 | Directory containing list of targets | Human editable
packages.bib                            | BibTex references for packages used  | Human editable
references.bib                          | BibTex references                    | Human editable
preregistration.rmd                     | Preregistered hypotheses             | Human editable
1_AEDB.CEA.baseline.Rmd                 | Preparing data, baseline table       | Human editable
3_1_bulkRNAseq.preparation.Rmd          | Preparing bulk RNAseq analyses       | Human editable
3_2_bulkRNAseq.exploration.Rmd          | Exploration RNAseq.                  | Human editable
4_scRNAseq.Rmd                          | Single-cell RNAseq analyses          | Human editable
6_Parsing_AE_molQTL.ipynb               | Parsing molQTL results               | Human editable
7_Parsing_GWASSumStats_MDD_vs_CAD.ipynb | Parsing GWAS SumStats MDD vs CAD     | Human editable
8_Extract_8_regions.ipynb               | Extract data for regions of interest | Human editable

<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

### Reproducibility

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to ensure transparency and reproducibility. The workflow is designed to meet the principles of Open Science throughout a research project. 

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles, read the preprint at https://osf.io/zcvbs/

#### WORCS: Advice for authors

* To get started with `worcs`, see the [setup vignette](https://cjvanlissa.github.io/worcs/articles/setup.html)
* For detailed information about the steps of the WORCS workflow, see the [workflow vignette](https://cjvanlissa.github.io/worcs/articles/workflow.html)

#### WORCS: Advice for readers

Please refer to the vignette on [reproducing a WORCS project]() for step by step advice.
<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->

# Acknowledgements
Dr. Sander W. van der Laan is funded through EU H2020 TO_AITION (grant number: 848146), EU HORIZON NextGen (grant number: 101136962), EU HORIZON MIRACLE (grant number: 101115381), and Health~Holland PPP Allowance ‘Getting the Perfect Image’.

We are thankful for the support of the Leducq Fondation ‘PlaqOmics’ and ‘AtheroGen’, and the Chan Zuckerberg Initiative ‘MetaPlaq’. The research for this contribution was made possible by the AI for Health working group of the [EWUU alliance](https://aiforhealth.ewuu.nl/). The collaborative project ‘Getting the Perfect Image’ was co-financed through use of PPP Allowance awarded by Health~Holland, Top Sector Life Sciences & Health, to stimulate public-private partnerships.

Plaque samples are derived from endarterectomies as part of the [Athero-Express Biobank Study](https://doi.org/10.1007/s10564-004-2304-6) which is an ongoing study in the UMC Utrecht. We would like to thank all the (former) employees involved in the Athero-Express Biobank Study of the Departments of Surgery of the St. Antonius Hospital Nieuwegein and University Medical Center Utrecht for their continuing work. Lastly, we would like to thank all participants of the Athero-Express Biobank Study; without you these kinds of studies would not be possible.

The framework was based on the [`WORCS` package](https://osf.io/zcvbs/).

## Disclosures
Dr. Sander W. van der Laan has received Roche funding for unrelated work.

<a href='https://uefconnect.uef.fi/en/group/miracle/'><img src='images/UEF_Miracle_Logo-07.png' align="center" height="75" /></a> <a href='https://www.to-aition.eu'><img src='images/to_aition.png' align="center" height="75" /></a> <a href='https://www.health-holland.com'><img src='images/logo_NL_HealthHollland_Wit-Oranje_RGB.png' align="center" height="35" /></a> <a href='https://www.nextgentools.eu'><img src='images/NextGen_1_Red.png' align="center" height="35" /></a> <a href='https://www.era-cvd.eu'><img src='images/ERA_CVD_Logo_CMYK.png' align="center" height="75" /></a> <a href=''><img src='images/leducq-logo-large.png' align="center" height="75" /></a> <a href='https://www.fondationleducq.org'><img src='images/leducq-logo-small.png' align="center" height="75" /></a> <a href='https://osf.io/zcvbs/'><img src='images/worcs_icon.png' align="center" height="75" /></a> <a href='https://doi.org/10.1007/s10564-004-2304-6'><img src='images/AE_Genomics_2010.png' align="center" height="100" /></a>

#### Changes log
    
    _Version:_      v1.4.2</br>
    _Last update:_  2025-07-15</br>
    _Written by:_   Sander W. van der Laan (s.w.vanderlaan-2[at]umcutrecht.nl).
    
    **MoSCoW To-Do List**
    The things we Must, Should, Could, and Would have given the time we have.
    _M_

    _S_

    _C_

    _W_

    **Changes log**
    * v1.4.2 Cleaning up. 
    * v1.4.1 Fixed some issues with writing the baseline tables. 
    * v1.4.0 Total re-organization. Updated to AEDB. Updated to RNAseq (deeper sequencing). Added baseline table for samples in eQTL analyses (b37 version).
    * v1.3.1 Fixed baseline table writing. Added additional saving options (raw, normalized, and log-transformed data) for the bulk RNAseq.
    * v1.3.0 Some script changes. Update to AEDB. Update to RNAseq (deeper sequencing). 
    * v1.2.2 Some organisational updates. 
    * v1.2.1 Fixed some references in the README.md. 
    * v1.2.0 There were some critical fixes in the way each notebook starts. Updates to the way the local system is set up and packages are loaded. Updates to functions. Updates to the organisation of the various notebooks. 
    * v1.1.0 Major update to WORCS system. 
    * v1.0.6 Small bug fixes. 
    * v1.0.5 Added png for overlap-figure.
    * v1.0.5 Removed obsolete references to objects.
    * v1.0.4 Fixed a mistake in the chr X sample-file creation. Now the order matches the chr X data.
    * v1.0.3 Fixed weight of files (limit of 10Mb per file for templates). Renamed entire repo.
    * v1.0.2 Added sex-specific .sample-files. Added GWASToolKit input-files.
    * v1.0.0 Initial version. Add 'plaque vulnerability index', Fixed baseline table, added codes, and results. Created sample-files. 

--------------

#### [Creative Commons BY-NC-ND 4.0](LICENSE)
##### Copyright (c) 1979-2025 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.
