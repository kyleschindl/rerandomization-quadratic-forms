Reproducibility for “A Unified Framework for Rerandomization using Quadratic Forms”

This repository contains all code, data processing steps, and documentation necessary to reproduce the results (figures, tables, and session information).


---
Requirements

- R (version 4.4.2 or later)
- renv for package management (see below)

---
Setup & Installation

1. Clone the repository
   git clone https://github.com/kyleschindl/rerandomization-quadratic-forms.git
   cd rerandomization-quadratic-forms

2. Restore R package environment
   Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')"
   Rscript -e "renv::restore(repos='https://cloud.r-project.org')"

3. Run the full analysis
   Rscript scripts/run_analysis.R

   This will:
   - Source helper functions (R/replication_functions.R)
   - Generate Figures 1–3 in output/figures/
   - Generate Table 2 and other CSV outputs in output/tables/

---
Directory Structure

├── README.md               
├── renv.lock               
├── .gitignore              
│
├── data/      
│
├── R/
│   └── replication_functions.R  
│
├── scripts/
│   ├── run_analysis.R      
│   ├── generate_figures.R  
│   └── generate_tables.R   
│
├── output/
│   ├── calculations/         
│   ├── figures/                       
│   └── tables/     
│
└── acc_form.pdf        

