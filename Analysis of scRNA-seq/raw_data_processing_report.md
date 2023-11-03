## Raw data processing 

To reproduce results of raw data processing I describe here all steps which I did in detail. Following data analysis you can find in `.ipynb` file.

Python packages used: `bioconda`, `pyroe`, `alevin-fry`,  `salmon`

1. According to the advice in [this course](https://ivanek.github.io/analysisOfGenomicsDataWithR/03_RNAseq_intro_html.html#1_Introduction_to_RNA_sequencing), I downloaded reference genome and gene annotation files
	* Genome sequence, primary assembly (GRCh38) (PRI) [Fasta](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz)  from [Gencode](https://www.gencodegenes.org/human/) , saved as `genome.fa`
	* Comprehensive gene annotation (CHR) [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) from [Gencode](https://www.gencodegenes.org/human/) , saved as `genes.gtf`
	* Barcode whitelist for 10x Chromium v2 chemistry `10xv2permit.txt`  from [here](https://umd.box.com/shared/static/jbs2wszgbj7k4ic2hass9ts6nhqkwq1p) 
    
<br>
    
2. Here it is needed to generate so-called the _splici_ index (spliced transcripts + introns) using `pyroe make-splici` as described in [Single-cell best practices](https://www.sc-best-practices.org/introduction/raw_data_processing.html#brief-discussion). 

```bash
# Activate conda
# Building the splici reference
pyroe make-splici fasta/genome.fa genes/genes.gtf 90 splici_rl90_ref

# Activate alevin-fry
conda activate af

# Index the reference
salmon index -t $(ls splici_rl90_ref/*\.fa) -i salmon_index -p 8
```

```bash
salmon alevin -i ref_data/salmon_index  -l ISR -1 seq_data/rna_reads_1.fastq.gz -2 seq_data/rna_reads_2.fastq.gz -o salmon_alevin -p 8 --chromium --sketch
```
the output is stored in `salmon-alevin/report_mapping.rtf`.

**Attention!** Here and further some filenames are changed for breavity. You can find the right filenames in `sensitive_info.txt`. 

<br>

3. Generate a permit list for cell barcode correction 

We have a white list for barcodes that allows us to search the cell barcodes in our data. By doing this, we can filter actual cells, which were probably in our sample. We also map to complementary DNA using the parameter `-d fw`. The alevin-fry output can be found in `alevin_fry_gpl/report_gpl.rtf` .

``` bash
# Cell barcode correction
alevin-fry generate-permit-list -u 10xv2permit.txt -d fw -i salmon_alevin -o alevin_fry_gpl
```

4.  Filtering 
Using the white list and the barcode mapping, we will put together the mapping of the sequencing reads against an index of the reference file (it generates `salmon_alevin/map.rad` file).
```bash 
alevin-fry collate -i alevin_fry_gpl -r salmon_alevin -t 8
```
the output of this command can be found in `alevin_fry_gpl/report_collate.rtf`. 

5. UMI resolution: allocate a molecular count to each gene in each cell.

`cr-like` creates a list of (gene, UMI, count) groups for each cell barcode. If a read matches multiple genes, it generates multiple groups. The groups are sorted based on gene, UMI, and count. UMIs matching one gene are assigned to that gene. UMIs matching multiple genes go to the gene with the highest count. If there's a tie, the corresponding reads are discarded.

At this stage, it is possible to use different strategies, for example, `cr-like-em`, instead of discarding UMIs with equal matches to multiple genes, it treats them as a group (a.k.a. equivalence class), using an expectation maximization algorithm to determine counts for each gene.

```bash 
## Usage: alevin-fry quant -r resolution -m txp_to_gene_mapping -i gpl_out_dir -o quant_out_dir -t num_threads

alevin-fry quant -r cr-like-em  -m $(ls ref_data/splici_rl90_ref/*3col.tsv) -i alevin_fry_gpl -o alevin_fry_quant -t 8
```
the report in the file `report_gpl.rtf`. 

One can map the RNA data with the same instructions as before for the other parts of data (only `fastq` file names have changed on `rna_reads_2_1.fastq.gz`, `rna_reads_2_2.fastq.gz`) or add modify the mapping step on 

```bash
salmon alevin -i ref_data/salmon_index  -l ISR -1 seq_data/rna_reads_1.fastq.gz, seq_data/rna_reads_2_1.fastq.gz  -2 seq_data/rna_reads_2.fastq.gz, seq_data/rna_reads_2_2.fastq.gz -o salmon_alevin -p 8 --chromium --sketch
```

## Antibody raw data processing

First of all, we need to index the feature barcodes. I received the barcode sequences of the antibody-derived tags (ADT) and the hash antibody tag oligos (HTO).


from `ref_data` folder

```bash
conda activate salmon
salmon index -t antibodies/adt_hto.tsv -i adt_hto_index --features -k7 
```
The output is in the `adt_hto_index` folder.
The command line output is stored in `adt_hto_index/report_adt_hto_mapping.rtf`

First, we map the ADT and HTO data. 

```bash
salmon alevin -l ISR -i ref_data/adt_hto_index -1 seq_data/antibodies/ant_reads_1.fastq.gz -2 seq_data/antibodies/ant_reads_2.fastq.gz --chromium -o adt_hto_mapping -p 8 --sketch
```

the mapping will just map each feature name to itself. 

```bash
awk '{print $1"\t"$1;}' adt_hto.tsv > t2g_adt_hto.tsv
```

Now, we quantify the ADT and HTO data
```bash
alevin-fry generate-permit-list -u 10xv2permit.txt -d fw -i adt_hto_mapping -o adt_hto_quant

alevin-fry collate -i adt_hto_quant -r adt_hto_mapping -t 8

alevin-fry quant -r cr-like-em -m ref_data/antibodies/t2g_adt_hto.tsv -i adt_hto_quant -o adt_hto_quant/crlike-em  -t 8 
```

We have mapped and quantified all of the data. Each of the relevant output directories `alevin_fry_quant_part_1/crlike` (for the first part of rna data), `alevin_fry_quant_part_2/crlike` (for the second part of rna data), `adt_hto_quant/crlike-em`. 

The command line output is stored in `adt_hto_quant/crlike-em/report_gpl_adt_hto.rtf`, `adt_hto_quant/crlike-em/report_collate_adt_hto.rtf`


**Attention!**  There was a problem with the detection of `Donor-4`. I suppose, that it was an edge effect: I added an empty row space after `Donor_4`, and the `af` assigned zero counts to the`' '` name, not for `Donor-4`. Now I need to delete the column with zero counts in the AnnData object, but at least I have counts for `Donor-4` now. 



---

#### Creation of a virtual environment 

```bash 

conda create --name bio_env_39 python=3.9
conda activate bio_env_39
conda install -c conda-forge scanpy python-igraph leidenalg
conda install pandas==1.5.3 --force-reinstall
conda install -c conda-forge anndata
conda install -c conda-forge pyroe

conda install -c conda-forge muon
conda install -c conda-forge mudata
conda install -c conda-forge scipy

#for doublet detection
conda install -c anaconda scikit-learn


conda install ipykernel
ipython kernel install --user --name=bio_env_39

```

