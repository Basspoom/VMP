<div align="center">


![VMP](VMP.jpg)

![biorxiv](https://img.shields.io/badge/bioRxiv-preprint-yellow.svg)
![HuggingFace Hub](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Models-red)
![Zenodo](https://img.shields.io/badge/zenodo-Datasets-green) 
![license](https://img.shields.io/badge/License-CC--BY--NC%204.0-blue.svg?labelColor=gray)

</div>

- [VMP](#VMP)
  - [Description](#description)
  - [Installation](#installation)
    - [Clone the source code](#clone-the-source-code)
    - [Create the conda environment](#create-the-conda-environment)
    - [Prepare required databases, VPAC models and the config file](#prepare-required-databases-VPAC-models-and-the-config-file)
  - [Quick Start](#quick-start)
    - [Automated Run](#automated-run)
    - [Quality Control](#quality-control)
    - [Assembly](#assembly)
    - [Viral Contig Identification](#viral-contig-identification)
    - [Clustering](#clustering)
    - [Binning](#binning)
  - [Citing this work](#citing-this-work)
  - [Acknowledgments](#acknowledgments)
  - [Contacts](#Contacts)




## Description

The **VMP (Virus Mining Pipeline)** is a modular and automated pipeline for mining viral genomes from raw microbiome data.
It integrates state-of-the-art tools for **quality control, assembly, viral contigs identification (our self-developed VPAC), clustering of viral genomes (vOTUs), and binning of microbial genomes (MAGs)**.
With its AI-empowered classifier VPAC, VMP provides both high accuracy and scalability across diverse microbiomes.


## Installation


### Clone the source code

```bash
git clone https://github.com/Basspoom/VMP.git
cd VMP
```


### Create the conda environment

The VMP uses three separate Conda environments to manage dependencies effectively. The `environment` files are located in the `VMP/environments` directory.

```bash
conda env create -f environments/VMP.yml -n VMP
conda env create -f environments/VPAC-single.yml -n VPAC-single
conda env create -f environments/VPAC-dual.yml -n VPAC-dual
```


### Prepare required databases, VPAC models and the config file

The VMP relies on several external databases and pre-trained models. Please make sure to download and configure them before the first run. The dependencies of each module are summarized below:

1. **Quality Control**: Host decontamination database (e.g., human reference genome hg38).
2. **Assembly**: No external databases required.
3. **Viral Contig Identification**: 
 - VIBRANT database: https://github.com/AnantharamanLab/VIBRANT/blob/master/
 - CheckV database: https://pypi.com.cn/project/checkv/
 - Kaiju database: https://github.com/bioinformatics-centre/kaiju/tree/master
 - geNomad database: https://github.com/apcamargo/genomad/
 - VPAC models: available via Zenodo https://zenodo.org/uploads/14033148  ***including single- and dual-path classifiers; Evo pre-trained models with path, model, and config***
4. **Clustering**: No external databases required.
5. **Binning**: No external databases required.

Users can run the script under the `VMP/database/` directory to automatically download and configure all required databases:
```bash
cd VMP/database
python prepare_database.py -db all
```

This script will fetch and set up all dependencies listed above.
- Estimated time: ~ xxx hours (depending on your network bandwidth).
- Disk space required: ~ xxx GB.
Typically, users only need to run this step once before the first use.


If some databases were already downloaded, users can specify only the missing ones. For example:
```bash
python prepare_database.py -db VirSorter2,CheckV,geNomad
```

After downloading, update the paths of databases in your VMP/config.yml file to ensure VMP can properly locate these databases and models, for example:
> ```yml
> paths:
>    VIBRANT_db: 'VMP/databases/VIBRANT_db'
>    CheckV_db: 'VMP/databases/CheckV_db'
>    Kaiju_db: '/your/path/to/Kaiju_db'
>    geNomad_db: '/your/path/to/geNomad_db'
>    VPAC_models: 'VMP/databases/VPAC_models'
> ```







## Quick Start

The GEM-PHI pipeline consists of three main scripts for feature calculation and one script for orchestration. All scripts are located in the `VMP/bin/` directory.

1. `calculate_node_features.py`: This script is the first step of the pipeline. It takes host and phage FASTA files as input and computes a variety of node-level features, including GC content, k-mer frequency, nucleotide transformer embeddings, and RBP embeddings for phages. This script must be run in the `GEM-PHI_nodes` Conda environment.

2. `calculate_edge_features.py`: The second step, this script calculates features for the edges connecting the nodes. It computes phage-phage, host-host, and phage-host edge features based on similarities and homologies. This script must be run in the `GEM-PHI_edges` Conda environment.

3. `final_inference.py`: This is the final step, which loads the pre-trained GEM-PHI model and the calculated node and edge features to predict phage-host interactions. The script outputs a list of predicted interactions and their confidence scores. This script must be run in the `GEM-PHI_inference` Conda environment.

For a complete and automated run of the entire pipeline, use the `run_all.py` script. It orchestrates the execution of the three sub-scripts in their respective Conda environments.


### Command Template

```bash
python ./bin/run_all.py \
  --config_path [path_to_config.yml] \
  --host_fasta [path_to_host_fasta] \
  --phage_fasta [path_to_phage_fasta] \
  --output_dir [output_directory] \
  --num_workers [number_of_cpu_cores] \
  --device [device_id]
```


### Run Example

```bash
python ./bin/run_all.py \
  --config_path ./config.yml \
  --host_fasta ./examples/hosts.fasta \
  --phage_fasta ./examples/phages.fasta \
  --output_dir ./outputs/example_run \
  --num_workers 16 \
  --device 0
```







## Citing this work
```bash
*@article {xxxx,
    author = {Bai, Zi-Peng and Zhao, Heng-Rui and Zhang, Ling-yu and Li, Bo-Rui and Chen, Yuxing and Qiong, Li and Zhou, Cong-Zhao},
    title = {VMP: an AI-empowered pipeline for mining viruses in microbiomes},
    elocation-id = {xxx},
    year = {2025},
    doi = {xxx},
    publisher = {xxx},
    URL = {xxx},
    eprint = {xxx},
    journal = {Microbiome}
}
```



## Acknowledgments

We acknowledge ***Hefei HiDimension Biotechnology Co., Ltd.*** for providing technical assistance and resources!!! 🥰🤗




## Contacts

We are honored to help you if you have any questions. Please feel free to open an issue or contact us directly. Hope our code helps and look forward to your citations.

[zcz@ustc.edu.cn] | [basspoom@mail.ustc.edu.cn] | [zhr123456@mail.ustc.edu.cn].
