# GSEA-SNP Workshop Tutorial

**Start here.** The original tutorial has been split so we can finish the essential `fgsea` exercise during class and leave the longer GenABEL/GWAA-to-GSEA project for independent exploration.

> **Acknowledgment:** Updated workshop adaptation of the **Neibergs Lab GSEA-SNP Tutorial**, originally written by **Tori Kelson and Allison Herrick**. 

## Which file should I use?

| File | Purpose | Required? |
|---|---|---|
| **[01_In_Class_Tutorial.md](01_In_Class_Tutorial.md)** | Built-in `fgsea` example; interpretation; plots; export filtered pathway results as CSV/XLSX; five-question checkpoint | **Yes** |
| **[02_Optional_Project.md](02_Optional_Project.md)** | GenABEL phenotype permutations; bovine GWAA-to-gene ranking; full GSEA; interpretation and troubleshooting | **No — optional reference** |

## What to do during class

1. Open **[01_In_Class_Tutorial.md](01_In_Class_Tutorial.md)** and follow the instructions in order.
2. Run the `fgsea` example and plot the enrichment results.
3. At Step 8, export the pathway results using **`NES > 2.5` and `padj < 0.05`**, then save the full results separately.
4. Answer the five questions under **🧠 Checkpoint**.
5. **STOP there.** Do **not** start GenABEL or the optional project during class.

## Files for homework

The required *filtered* result export is created in:

```text
~/workshop/gsea/fgsea_example/
```

Expected files:

```text
fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.csv
fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.xlsx
```

The Excel file requires the R package `writexl`; **the CSV contains the same results and can be used if `writexl` is unavailable**. Students also create an unfiltered `fgsea_example_results.csv` to support broader comparisons and checkpoint questions.

Each row of the filtered file is a pathway with **positive NES above 2.5** and **adjusted p-value below 0.05**, together with `all_genes`, `leadingEdge`, and `leading_edge_count` columns. The filter is intentionally positive-only; negative NES pathways will not appear.

## Later, if you are curious

Go to **[02_Optional_Project.md](02_Optional_Project.md)**. Its major sections are collapsed; expand only those you want to study. The optional project uses GWAA results and shared files and is **not needed** to complete this week's example-based homework files.

---

*Workshop topic: Gene Set Enrichment Analysis (GSEA-SNP) following bovine GWAA.*
