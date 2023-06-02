# Variant catalogue Pipeline


## Introduction

**The variant catalogue Pipeline** is a workflow designed to generate variant catalogues, a list of variants and their frequencies in a population, from whole genome sequences.

the variant catalogue pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It takes as input Whole Genome Sequence (WGS) data and outputs multiple vcf files including the variant allele frequencies in the cohort and some basic annotation.

The variant catalogue pipeline includes detection of Single Nucleotide Variants (SNV), small insertions and deletions (indels), Mitochondrial variants, Structural Variants (SV), Mobile Element Insertions (MEI), and Short Tandem Repeats (STR). The output variant catalogue can be generated for GRCh37 and/or GRCh38 human reference genomes.

<p align="center">
    <img title="The variant catalogue Workflow" src="https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/supplementary_information/Variant_catalogue_figure.png" width=50%>
</p>
<p align="center">
Figure : Overview of the variant catalogue pipeline
</p>


## Pipeline description

The variant catalogue is composed of four sub-workflows represented by the grey boxes in the figure. The structure allows users to run the pipeline as a whole or choose to run individual sub-workflow(s) of interest. Each sub-workflow is composed of modules, which call upon open-access genomic software.

A more detailed description of the pipeline will be available soon.

Detailed representation of the variant catalogue pipeline : [Variant_catalogue_supp_figure.pdf](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/supplementary_information/Variant_catalogue_supp_figure.pdf)

## Pipeline availability

The variant catalogue pipeline is implemented in the NextFlow framework and relies only on open-access tools, therefore, any user with sufficient compute capacity should be able to use this pipeline. Users who want to use this pipeline on their local servers will have to install the necessary software on their instance.

All the software required to run the variant catalogue pipeline are open-source and the link to the installation guidelines are available in [supplementary_information/software_information.md](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/software_information.md).

All the other resources necessary to run the pipeline (Reference genomes, annotation plugins, etc) are also publicly available and information related to them are available for [GRCh37](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/GRCh37_specific_files.md) in supplementary_information/GRCh37_specific_files.md, for [GRCh38](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/GRCh38_specific_files.md) in supplementary_information/GRCh38_specific_files.md and for [the mitochondrial genome](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/Mitochondrial_references.md) in supplementary_information/Mitochondrial_references.md.

## Step-by-step instructions to run the pipeline

**1. Download and install the softwares necessary to run the pipeline.**

All the software required to run the variant catalogue pipeline are open-source and the link to the installation guidelines are available in [supplementary_information/software_information.md](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/software_information.md).

**2. Download all the resources necessary to run the pipeline** (Reference genomes, annotation plugins, etc).

For [GRCh37](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/GRCh37_specific_files.md) in supplementary_information/GRCh37_specific_files.md

For [GRCh38](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/GRCh38_specific_files.md) in supplementary_information/GRCh38_specific_files.md

For [the mitochondrial genome](https://github.com/wassermanlab/CAFE_pipeline/blob/main/supplementary_information/Mitochondrial_references.md) in supplementary_information/Mitochondrial_references.md.

**3. (Optional) Download a set of publicly available genomes.
Or use your own set of genomes.**

The samples that were processed in the test dataset are available here : [batch_1](https://github.com/wassermanlab/CAFE_pipeline/blob/main/test_case/80_samples_information) contained 80 samples and [batch_2](https://github.com/wassermanlab/CAFE_pipeline/blob/main/test_case/20_samples_information) contained 20 samples.

**4. Clone this repository.**

```
git clone wassermanlab/Variant_catalogue_pipeline
```

**5. Update the config file and the bash file.**

Upade the path in the nextflow.config and the Variant_Catalogue_v1.sh files to point to the software and files as per your organization.

The pipeline is set up to run using SLURM scheduler, however, Neftflow is adaptable to other system. If another system is used, the nextflow.config file should be updated.

For more information regarding nextflow config file, visit [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html)

If you want to run the pipeline with GRCh38 as a reference, replace GRCh38 by GRCh37 in the Variant_Catalogue_v1.sh file

**6. Launch the pipeline.**

```
bash Variant_Catalogue_v1.sh
```


## Future of the pipeline

There is discussions with the [NF_core](https://nf-co.re) team to move the variant catalogue pipeline to NF-core and improve and maintain this pipeline with the help of the community. Feel free to join in the effort!


## Pipeline test on 100 genomes

In order to test the variant catalogue pipeline, 100 samples from the IGSR (International Genome Sample Resource) were processed. A more precise description of the method and results will be available soon. 

The samples were processed in two batches, [batch_1](https://github.com/wassermanlab/CAFE_pipeline/blob/main/test_case/80_samples_information) contained 80 samples and [batch_2](https://github.com/wassermanlab/CAFE_pipeline/blob/main/test_case/20_samples_information) contained 20 samples.
Output files generated by Nextflow (report, timeline, etc) are available in the test_case folder and vcf files containing annotated variant frequencies are also available in that folder. Intermediate files were not loaded into GitHub.

## Pipeline contributors

[Solenne Correard](https://github.com/scorreard) : Developed the SNV and MT sub-workflows, contributed to hail script and other downstream modules, lead the testing of the pipeline.

[Mohammed OE Abdallah](https://github.com/melsiddieg) : Developed the SV sub-workflow, contributed to hail setup and script and to the migration to Nextflow DSL2.

[Brittany Hewitson](https://github.com/brittanyhewitson) : Reviewed the modules and sub-workflows and contributed to the final outputs layout.

[Phillip Richmond](https://github.com/Phillip-a-richmond) : Contributed to testing and configuration.

