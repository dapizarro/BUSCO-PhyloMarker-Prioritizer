# PhyloMarker Prioritizer

**BUSCO-derived prioritization of phylogenetic markers from genome-scale datasets.**

This repository implements a reproducible pipeline to identify, rank, and select informative single-copy orthologous genes for phylogenomic analyses. It was designed for lichenized fungi and Pezizomycotina-scale datasets, but the workflow is general and can be applied to any taxonomic group with BUSCO outputs.

## Scientific rationale

Genome-scale phylogenetics often starts from hundreds or thousands of orthologous loci. BUSCO is useful because it identifies broadly conserved, near-universal single-copy genes, but **single-copy orthology does not automatically mean high phylogenetic utility**. BUSCO genes differ in evolutionary rate, taxonomic occupancy, alignment quality, proportion of missing data, and number of phylogenetically informative sites.

This pipeline prioritizes loci using three complementary dimensions:

1. **Taxonomic representation**: genes present across more taxa are less affected by missing data and provide a more stable backbone.
2. **Phylogenetic informativeness**: genes with more variable and parsimony-informative sites are more likely to resolve relationships.
3. **Alignment quality**: genes with excessive gaps or ambiguous regions can introduce noise and artefacts.

The goal is not to find one universally best gene, but to create **ranked marker panels** such as Top25, Top50, Top100, Top164, Top200, or a data-driven optimum panel. These panels can then be tested by concatenation, gene-tree/species-tree inference, concordance analysis, or targeted bait design.

## Conceptual workflow

```text
BUSCO outputs per genome
        │
        ├── BUSCO summary table
        │
        ├── Extract single-copy BUSCO sequences
        │        one FASTA per BUSCO gene
        │
        ├── MAFFT protein alignment
        │
        ├── trimAl conservative trimming
        │
        ├── Per-gene metrics
        │        occupancy, gap fraction, variable sites,
        │        parsimony-informative sites, entropy
        │
        ├── PCA-based ranking
        │        PC1 = synthetic marker quality axis
        │
        ├── Saturation / breakpoint analysis
        │        optimal panel size
        │
        └── Export ranked marker panels and diagnostic plots
```

## Input data

The expected input is a directory containing one BUSCO run per genome:

```text
data/busco_runs/
├── busco_Taxon_1/
│   └── run_ascomycota_odb10/
│       └── busco_sequences/
│           └── single_copy_busco_sequences/
│               ├── 1001at4890.faa
│               └── ...
├── busco_Taxon_2/
└── busco_Taxon_3/
```

Only complete single-copy BUSCO sequences are used by default. Multi-copy and fragmented BUSCOs are deliberately excluded because the first aim is to build a clean orthologous marker framework.

## Installation

Using conda/mamba:

```bash
mamba env create -f environment.yml
mamba activate phylo-marker-prioritizer
```

Or using pip for the Python dependencies, assuming `mafft` and `trimal` are already installed:

```bash
pip install -r requirements.txt
```

## Quick start

From the repository root:

```bash
bash run_pipeline.sh data/busco_runs results
```

For a cluster or workstation, set the number of MAFFT threads per gene:

```bash
THREADS=4 bash run_pipeline.sh data/busco_runs results
```

## Main outputs

```text
results/
├── busco_summary_all.tsv
├── busco_per_gene/
│   └── <BUSCO_ID>.faa
├── alignment_per_gene_trimmed/
│   └── <BUSCO_ID>.trim.faa
├── per_gene_metrics.tsv
├── per_gene_metrics_PCA.tsv
├── optimal_panel_optimal_panel_size.tsv
├── optimal_panel_cumulative_curves.tsv
├── optimal_panel_optimal_panel_curve.png
├── plots/
│   ├── corr_heatmap.png
│   ├── hist_*.png
│   ├── scatter_PC1_vs_*.png
│   └── tau_by_metric.png
└── panels/
    ├── top25_genes.txt
    ├── top50_genes.txt
    ├── top100_genes.txt
    └── ...
```

The most important output is `per_gene_metrics_PCA.tsv`, which contains per-gene metrics, PC1/PC2 scores, and the final ranking.

## Metrics calculated per gene

| Metric | Meaning | Desired direction |
|---|---|---|
| `n_taxa` | Number of taxa represented in the alignment | High |
| `aln_len` | Trimmed alignment length | Context-dependent |
| `occupancy` | Non-missing cells / total cells | High |
| `gap_fraction` | Gap characters / total cells | Low |
| `gap_cols_fraction` | Fraction of columns containing at least one gap | Low |
| `variable_sites` | Number of variable alignment columns | High |
| `PIS` | Parsimony-informative sites | High |
| `VS_per_len` | Variable sites divided by alignment length | High |
| `PIS_per_len` | PIS divided by alignment length | High |
| `PIS_over_VS` | Fraction of variable sites that are parsimony-informative | High |
| `mean_entropy` | Mean Shannon entropy across informative columns | Moderate to high |

## Ranking strategy

The default ranking uses PCA on standardized metrics:

```text
occupancy
VS_per_len
PIS_per_len
PIS_over_VS
mean_entropy
```

PC1 is interpreted as a synthetic axis of marker quality. Genes with higher PC1 values are ranked higher. This avoids relying on one arbitrary scalar score and allows several properties of marker quality to contribute jointly.

The PCA ranking is useful because marker quality is multidimensional: a gene may have many informative sites but poor occupancy, or excellent occupancy but too little variation. PCA provides an objective, reproducible compromise among these dimensions.

## Optimal panel size

The script `find_optimal_panel_size.py` evaluates cumulative median behaviour along the ranked gene list. For each TopN panel, it calculates cumulative metric curves and estimates a two-segment breakpoint, interpreted as a saturation point. The final recommended panel size is the median breakpoint across metrics.

This gives a data-driven answer to the question:

> How many genes are enough before adding more loci gives diminishing returns?

The optimal panel is not a replacement for biological validation. It should be compared against fixed-size panels and full datasets using downstream phylogenetic analyses.

## Recommended downstream validation

After generating ranked panels, test them with:

1. Concatenated ML trees with partitioning/model selection, e.g. IQ-TREE2.
2. Gene-tree inference for each selected locus.
3. Species-tree inference, e.g. ASTRAL.
4. RF distances among panels.
5. Concordance factors, quartet support, and conflict analysis.
6. Focused testing in difficult clades or recently diversified lineages.

## Example commands for downstream phylogenomics

```bash
# Example: concatenate selected alignments with AMAS.py after copying TopN loci
AMAS.py concat -i selected_alignments/*.trim.faa \
  -f fasta -d aa \
  -t results/panels/top100_supermatrix.fas \
  -p results/panels/top100_partitions.txt

# IQ-TREE2 example
iqtree2 -s results/panels/top100_supermatrix.fas \
  -p results/panels/top100_partitions.txt \
  -m MFP+MERGE -B 1000 -alrt 1000 -T AUTO
```

## Citation / conference abstract

If using this workflow, cite the associated conference contribution:

> Pizarro D., Vaiana A., Calzoni D., Singh G., Lumbsch T. & Divakar P.K. Which genes matter? Prioritizing phylogenetic signal in lichenized fungi.

## Repository structure

```text
phylo_marker_prioritizer/
├── README.md
├── environment.yml
├── requirements.txt
├── run_pipeline.sh
├── config/
│   └── config.yaml
├── scripts/
│   ├── make_busco_summary_table.py
│   ├── extract_singlecopy_busco_per_gene.py
│   ├── align_and_trim.py
│   ├── compute_alignment_metrics.py
│   ├── rank_genes_pca.py
│   ├── find_optimal_panel_size.py
│   ├── plot_metrics.py
│   └── export_panels.py
├── docs/
├── tests/
└── .github/workflows/
```

## Notes and limitations

- BUSCO markers are conserved by design; they are ideal for backbone and comparative marker selection, but may not always resolve very shallow radiations.
- Protein alignments are recommended for deep phylogenomic marker discovery; nucleotide/CDS recovery can be performed later for shallower analyses.
- The ranking is dataset-dependent. A Top100 panel optimized for Pezizomycotina should be re-evaluated before being transferred to another lineage.
- PCA ranking identifies a quantitative compromise, not a biological truth. Always validate the resulting panels phylogenetically.

## License

MIT.
