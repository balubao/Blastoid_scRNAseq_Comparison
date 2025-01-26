# Blastoid_scRNAseq_Comparison

This repository contains the scripts and workflows used in the analysis presented in Balubaid et al., 2024, *iScience*. The study evaluates single-cell RNA sequencing (scRNAseq) data from blastoid models.

## Overview
The analysis was conducted through the following steps:

1. **Generating Count Matrices:**
   - Processing FASTQ files to generate count matrices using standard pipelines.

2. **Preprocessing Data:**
   - Quality control (QC) and filtering of count matrices to generate high-quality data for downstream analysis.

3. **Cell Type Annotation:**
   - Annotating cell types using a consensus annotation scheme across datasets.

4. **Dataset Integration:**
   - Integrating multiple datasets to enable comparative analysis.

5. **Data Analysis:**
   - **Cluster Distribution Analysis:** Identifying and characterizing cell clusters.
   - **Lineage Module Score Analysis:** Assessing lineage-specific gene expression patterns.

6. **Founding Cell Line Analysis:**
   - Sequencing and analyzing founding cell lines to understand their contribution to blastoid models.

## Repository Structure
- `scripts/`: Contains all analysis scripts.
- `data/`: Placeholder for input data files (not included due to size and privacy restrictions).
- `results/`: Contains processed data, figures, and other output files.

## Dependencies
This analysis was performed using the following tools and libraries:

- **R (v4.4.1):** R Core Team (2017). R: A language and environment for statistical computing. [R Project](http://www.r-project.org/); RRID:SCR_001905
- **Seurat (v5.0.1):** Hao et al. [Seurat](https://satijalab.org/seurat/)
- **SCINA (v1.2.0):** Zhang et al. [SCINA](https://github.com/jcao89757/SCINA)
- **SCType:** Ianevski et al. [SCType](https://github.com/IanevskiAleksandr/sc-type)
- **SingleR (v2.6.0):** Aran et al. [SingleR](https://www.bioconductor.org/packages/release/bioc/html/SingleR.html)
- **Seurat-Disk:** [GitHub](https://github.com/mojaveazure/seurat-disk)
- **Seurat-Wrapper:** [GitHub](https://github.com/satijalab/seurat-wrappers)
- **Harmony (v1.2.1):** Korsunsky et al. [Harmony](https://github.com/immunogenomics/harmony)
- **Stats:** R Core Team (2017). [Stats](https://www.r-project.org/)
- **philentropy (v0.8.0):** Drost et al. [philentropy](https://cran.r-project.org/web/packages/philentropy/)
- **Pheatmap (v1.0.12):** Kolde [Pheatmap](https://cran.r-project.org/web/packages/pheatmap/)
- **AUCell (v1.26.0):** Aibar et al. [AUCell](https://www.bioconductor.org/packages/release/bioc/html/AUCell.html)

For a complete list of dependencies, refer to `environment.yml`.

## Reproducibility
To reproduce the analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/<username>/Blastoid_scRNAseq_Comparison.git
   cd Blastoid_scRNAseq_Comparison
   ```

3. Run the scripts in the `scripts/` directory following the order specified in the documentation.

## Citation
If you use this repository, please cite the following:

Balubaid et al., 2024, *iScience* [link](https://www.cell.com/iscience/fulltext/S2589-0042(24)02347-2).

## Contact
For questions or feedback, please get in touch with Ali Balubaid at ali.balubaid(at)kaust.edu.sa.


