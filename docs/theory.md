# Theoretical background

This workflow treats marker selection as a multi-objective prioritization problem. A useful phylogenetic marker should be orthologous, broadly sampled, alignable, and sufficiently variable to distinguish evolutionary histories. BUSCO provides a practical orthology framework, but the final marker panel should still be evaluated quantitatively because BUSCO loci are heterogeneous in missing data, evolutionary rate, and alignment quality.

The implemented ranking combines occupancy, variable-site density, parsimony-informative-site density, the ratio of informative to variable sites, and entropy. These metrics are standardized and summarized by PCA. The first principal component is used as the default ranking axis because it captures the dominant joint gradient of marker quality in the dataset.

The optimal panel size is estimated from cumulative saturation curves. As ranked genes are added one by one, the median value of each metric is recalculated. A two-segment breakpoint model identifies where the improvement curve changes slope. The consensus of these breakpoints provides an empirical panel size for downstream testing.
