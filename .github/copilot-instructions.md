# Copilot Instructions for htlv_paper_2025

## Project Overview
This repository contains data analysis, scripts, and manuscript materials for the 2025 HTLV research paper. The project is primarily R-based, with some Python and shell scripts, and includes Quarto/Markdown documents for reporting and manuscript generation.

## Directory Structure
- `2503_bash_scripts/`: Shell scripts for file renaming and merging.
- `250212_clonality_table1/`: R scripts, Quarto files, and data for clonality analysis.
- `250304_HLAI/`: HLA analysis (R, Python, data, and results).
- `2506_mutations/`: Mutation analysis scripts, images, and raw data.
- `2507_PVL/`: PVL analysis scripts and Quarto files.
- `images/`, `output/`: Figures, tables, and supplementary outputs.
- `Manuscript.qmd`, `Manuscript.html`, `250508_Manuscript.docx`: Main manuscript files.

## Coding Conventions
- Use R (tidyverse style) for data wrangling, analysis, and visualization.
- Use Quarto (`.qmd`) for reproducible reports and manuscript generation.
- Use `maftools`, `ggplot2`, `pheatmap`, `flextable`, `officer`, and `kableExtra` for analysis and reporting.
- Python scripts (e.g., `hla_hed.py`) are used for specialized analyses.
- Shell scripts are for file management and preprocessing.

## AI Agent Guidance
- **Prioritize R and Quarto workflows** for mutation, HLA, and PVL analyses.
- **Do not overwrite raw data files** in any subdirectory.
- **When generating new scripts or reports:**
  - Place R scripts in the relevant analysis subfolder (e.g., `2506_mutations/`).
  - Place Quarto/Markdown reports in the same folder as the analysis or in the project root if general.
  - Place figures and tables in `output/` or the relevant analysis subfolder.
- **For statistical modeling:**
  - Use Poisson, negative binomial, beta regression, robust regression, or logistic regression as appropriate.
  - Document model assumptions and results clearly in Quarto reports.
- **For data wrangling:**
  - Use tidyverse idioms (e.g., `dplyr`, `tidyr`).
  - Clean and merge data frames before analysis.
- **For reporting:**
  - Use `flextable`, `officer`, or `kableExtra` for tables.
  - Export publication-ready figures (PNG, PDF) to `output/`.
- **For code suggestions:**
  - Follow existing file and folder naming conventions.
  - Use clear, descriptive variable and function names.
  - Add comments to explain non-obvious code.

## Special Instructions
- Do not modify manuscript files (`Manuscript.qmd`, `.docx`, `.html`) unless explicitly requested.
- Do not commit large raw data files or outputs to version control.
- When in doubt, ask for clarification before making major changes.

---
_Last updated: 2024-06_
