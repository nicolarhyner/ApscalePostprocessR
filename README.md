# ApscalePostprocessR  
### R scripts to refine and filter Apscale ESV tables

This repository provides an R script to **create, filter, and refine ESV tables** generated from the combined outputs of **APSCALE** and **apscale_blast / BOLDigger** metabarcoding pipelines.

> ⚠️ This is **not** an official APSCALE extension.  
> It is a helper script to post‑process Apscale outputs.

The script is **not fully automated and requires manual revision** before use.  
In particular, you may need to adjust column names to match your Apscale output files, create a vecor with your negative controlls and define the variables.

The current version was tested with:  
- **Apscale v2.1.1**  
- **apscale_blast v1.3.2**  
- **BOLDigger v3 (2.2.0)**  

All steps in the script are annotated, although some may still require additional clarification.

---

## Features
- Merge the ESV table with taxonomic assignments  
- Minimum abundance filter ESV- and sample wise
- Perform blank corrections using negative controls  
- *(Optional)* Collapse ESVs assigned to the same species  

After each step, an intermediate ESV table is exported.

---

## Required folder structure

The script is intended to run inside an **RStudio Project** with the following folder structure:
<Project_Name>
/Intermediate   # intermediate ESV tables
/Output         # final processed tables
/Scripts        # store the .R script here

These folders are created automatically by the script.

---

## ⚠️ To‑Do / Missing features

- Improve annotation and documentation  
- Add a Markdown report summarizing filters applied and abundances removed  
  (currently, these stats are printed only in the console)
- If you encounter issues / bugs, please let me know

---
