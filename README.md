# Script to convert GFF3 file to maizeGDB specifcations

MaizeGDB will assign the annotation prefix as well as provide you the genome/species accession ID. Once you have that you can run the script to get your GFF3 as per maizeGDBs requirements.


## Installation

The script itself does not require any extra packages installed, but since we have various formats of GFF3 files out there, we will need to standardize it before you can run this script - which requires you to install [`AGAT`](https://github.com/NBISweden/AGAT). Below is the recommended way to install this:


```bash
module load miniconda3
conda create -y -n agat
conda activate agat
conda install -c bioconda agat
```

Get the script:

```bash
git clone git@github.com:HuffordLab/MaizeGDB_gff3_format.git
chmod +x maizegdb_gff3_formatter.py
```

## Running the script

You can run the script as follows:

### Usage:

```
maizegdb_gff3_formatter.py <renamed_agat_formatted.gff3> <canonical_transcript_ids.txt> 
```


### Input:

`renamed_agat_formatted.gff3` : The gff3 file sanitized using the `agat_sp_manage_IDs.pl` script. Typically, you should request obtain gene id prefix from maizeGDB, and then run this as follows:

```
agat_sp_manage_IDs.pl --gff input.gff3 --prefix Ab00001aa --tair  --output prefinal.gff3
```

`canonical_transcript_ids.txt` : list of transcript ids that are considered as primary transcript. You can run the TRaCE program to determine the canonical transcript and create a list of mRNA ids (one per line). The number should be equal to the gene count in GFF3


### Example:

```
maizegdb_gff3_formatter.py renamed_agat_formatted.gff3 canonical_transcript_ids.txt > maizeGDB_specifications.gff3 
```



