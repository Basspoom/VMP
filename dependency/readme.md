# Prepare required databases, tools and models

The VMP relies on several external databases, several tools, and pre-trained models. Please make sure to download and configure them before the first run. The dependencies of each module are summarized below:

1. **Quality Control**: 
 - Host decontamination database (e.g., human reference genome hg38: http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz).
 - GenomeScope tool: https://github.com/schatzlab/genomescope/
2. **Assembly**: No external databases required.
3. **Viral Contig Identification**: 
 - DeepVirFinder tool: https://github.com/jessieren/DeepVirFinder/ 
 - VirSoeter2 database: https://osf.io/v46sc/download/ 
 - VIBRANT database: https://github.com/AnantharamanLab/VIBRANT/blob/master/
 - CheckV database: https://pypi.com.cn/project/checkv/
 - Kaiju database: https://github.com/bioinformatics-centre/kaiju/tree/master
 - geNomad database: https://github.com/apcamargo/genomad/
 - VPAC models: available via Zenodo https://zenodo.org/records/17008739  
 - Evo model: available via Zenodo https://zenodo.org/uploads/14033148 
4. **Clustering**: No external databases required.
5. **Binning**: No external databases required.

**Users can run the script under the `VMP/dependency/` directory to automatically download and configure all required files:**
```bash
cd VMP/dependency
python prepare_database.py -db all
```

This script will fetch and set up all dependencies listed above.
- Estimated time: 5~10 hours (depending on your network bandwidth).
- Disk space required: ~ 130 GB.
- Typically, users only need to run this step once before the first use.


**Optional:**   Enable Faster Downloading with aria2c
If available, the script can use `aria2c` for parallel accelerated downloading (no root required):
```bash
conda install -c conda-forge aria2
python prepare_database.py -db all -aria2
```

**Optional:**   Remove compressed packages after installation
Add the `-clean` flag:
```bash
python prepare_database.py -db all -aria2 -clean
```

If some databases were already downloaded, users can specify only the missing ones. For example:
```bash
python prepare_database.py -db VirSorter2,CheckV,geNomad --aria2
```