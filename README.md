# Positive selection analysis using Hyphy

### Install the tools needed
```
conda install -c bioconda minimap2
conda install -c bioconda samtools
conda install -c bioconda bcftools
```

### 1. Extraction of sequences of a specific gene
- Downloading the genome sequence `fasta` and annotation `gff3` of the reference `NC_045512.2`

  ```
  wget -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta" > NC_045512.2.fasta

  wget -O - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz" | gunzip > NC_045512.2.gff
  ```

- Extract the `Spike gene` sequence from the reference genome
  ```
  samtools faidx NC_045512.2.fasta NC_045512.2:21563-25384 > NC_045512.2.spike.fasta
  ```

- Map the sequences to the Spike gene sequence then sort and index
  ```
  minimap2 -ax asm5 NC_045512.2.spike.fasta alpha_550.fasta | samtools sort -o alpha.spike.sorted.bam
  
  samtools index alpha.spike.sorted.bam
  ```

- Convert the `bam` file which contains the sequences aligned to the `Spike gene` into a `fasta` format
  ```
  samtools fasta alpha.spike.sorted.bam > alpha.spike.fasta
  ```

### 2. Codon-aware Multiple Sequence Alignment
- Correct for frame-shift mutations
```
hyphy pre-msa.bf --reference NC_045512.2.spike.fasta --input alpha.spike.fasta
```

- Generate multiple sequence alignment of the spike protein sequence
```
muscle -align alpha.spike.fasta_protein.fas -output alpha.spike.fasta_protein.msa
```

- Generate corrected aligned nucleotide sequence of the spike gene 
```
hyphy post-msa.bf --protein-msa alpha.spike.fasta_protein.msa --nucleotide-sequences alpha.spike.fasta_nuc.fas --output alpha.fin.msa
```

### 3. Phylogenetic tree construction: Maximum Likelihood method
- Construct tree using maximum likelihood with `GTR` substitution model, `+I+G` invariable site plus discrete Gamma model, and `1000 bootstrapping`
```
iqtree2 -s alpha.fin.msa -m GTR+I+G -T AUTO -B 1000
```

### 4. Selection analysis
- Test for lineage-specific evolution using the branch-site method `aBS-REL` or `Adaptive Branch-Site Random Effects Likelihood`
```
hyphy absrel --alignment alpha.fin.msa --tree alpha.fin.msa.treefile CPU=10
```
