### 1. Concatenates input files and generate per read stats.

The [fastcat](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Filter out host reference reads

If a host_reference fasta file is provided, [Minimap2](https://github.com/lh3/minimap2) is used to align all reads to the host_reference, and any aligned reads are filtered out.

### 3. Trim reads

The reads are then trimmed at the ends using [SeqKit](https://bioinf.shenwei.me/seqkit/) with the provided trim length parameter, which has a default of 150bps. This default was chosen as we found it to produce optimum results and we don’t recommend changing this. At this stage SeqKit is also used to filter out reads that are greater than 1.2 x the approximate size or less than 100bps.

### 4. Subsample reads

The sequences are then subsampled using [Rasusa](https://github.com/mbhall88/rasusa). The subsampling will take the expected coverage parameter in to account; as we will be repeating the assembly 3 times, we subsample to a target of approximately 3x the required coverage. However, this is just a target and if there is not enough data, Rasusa is still able to create the 3 subsamples. The approximate size parameter is also used by Rasusa to work out the target number of bases and therefore number of reads required for each of the subsamples.

### 5. Create 3 subsamples

+[Trycycler](https://github.com/rrwick/Trycycler) is used to create 3 subsamples as we will be creating three assemblies and finding the consensus between all three. This consensus generation will be handled by Ttrycycler.

### 6. Assembly

We perform the assembly for each of the 3 subsamples separately using the assembly tool [Flye](https://github.com/fenderglass/Flye).

### 7. De-concatenate

If there are concatemers in the assembly, these are found using minimap2 and de-concatenated using a custom Python script. If the assembly is already roughly the expected approximate size, this de-concatenate step will be skipped.

### 8. Reconcile and polish

Trycycler is used to reconcile the subsampled assemblies into one final assembly. This is then polished with [Medaka](https://github.com/nanoporetech/medaka). A per-base quality score for the assembly is output by Medaka in a FASTQ file. This is used for creating the mean assembly quality you will find in the report.

### 8. Insert location and QC

SeqKit is used to locate inserts using the primers supplied to the primers parameter.

A multiple sequence alignment (MSA) will be done using [Pyspoa](https://github.com/nanoporetech/pyspoa). This will be presented in the report to help users compare inserts across samples in a multi-sample run. If an insert reference FASTA file is provided, this will also be included in the MSA.

If a reference insert FASTA sequence is provided, [BCFtools](https://samtools.github.io/bcftools/bcftools.html) is used to find variants between it and the final insert assembly, and are reported in BCF file per sample.

### 10. Annotate

The assembly is annotated by pLannotate](https://github.com/barricklab/pLannotate) to show any features that are present. The default database is used, which contains entries from [FPbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot), [Rfam](https://rfam.org/) and [SnapGene](https://www.snapgene.com/). Descriptions, percentage match and length of the match are also provided.

### 11. Self alignment

For each sample a self alignment will be done using [Last](https://gitlab.com/mcfrith/last) and the output will be presented as a dotplot. This can help identify any repetitive regions in your final assembly.