# 01 — In-Class Tutorial: Introduction to GSEA with `fgsea`

This tutorial continues directly from the **GWAA workshop tutorial**. You will first use the built-in `fgsea` example data to learn the mechanics of gene set enrichment analysis, then connect those concepts to the bovine GWAA dataset used in the workshop.

> **Acknowledgment:** This tutorial is an updated workshop adaptation of the **Neibergs Lab GSEA-SNP Tutorial**, originally written by **Tori Kelson and Allison Herrick**. Their tutorial established the two-part GenABEL/GSEA-SNP workflow, file-preparation logic, permutation strategy, and pathway-analysis framework used as the foundation for this updated version.

<sub>*For shorter commands, it is better to type the command manually rather than copying it directly. This helps you become more comfortable working in Bash and R.*</sub>

---

## Today’s learning goals (required)

In class, we will use the built-in `fgsea` example data to:

1. Inspect ranked genes and pathway gene sets.
2. Run GSEA and interpret `ES`, `NES`, `pval`, `padj`, and `leadingEdge`.
3. Make enrichment plots and examine the results.
4. **Export the filtered pathway results needed for homework** to CSV and Excel (if the Excel-writing package is available).
5. **Stop after the five checkpoint questions.**

> **Not required in class:** GenABEL, phenotype permutations, SNP-to-gene mapping, and running GSEA on the bovine GWAA results are retained in [02_Optional_Project.md](02_Optional_Project.md) for students who want to explore later. You do **not** need to complete these sections for Tuesday’s class or to create the example-results homework files.

---

# 💻 Server Login

You already completed the computer setup during the GWAA tutorial, so we will **not repeat the VPN, WSL, XQuartz, PuTTY, or FileZilla installation instructions here**.

> ⚠️ **Reminder:** If you are not connected to the WSU network, connect to the **WSU VPN** before logging in.

Use your assigned workshop username:

```bash
ssh -Y studentXX@10.104.58.24
```

For example:

```bash
ssh -Y student01@10.104.58.24
```

After logging in:

```bash
hostname
whoami
pwd
```

You should be on the workshop server and in your Linux home directory.

Now return to the same workshop directory used for the GWAA tutorial:

```bash
cd ~/workshop
pwd
ls
```

You should still see files generated during the GWAA tutorial, such as:

```text
srd_qc_hwe.bed
srd_qc_hwe.bim
srd_qc_hwe.fam

gwas_unadjusted_ADD.assoc.logistic
gwas_pc1_pc2_ADD.assoc.logistic
```

> 💡 The **in-class `fgsea` example does not use these GWAA association files**. They are only needed for the [optional GWAA → GSEA project](02_Optional_Project.md).

Create a separate directory for the GSEA work:

```bash
mkdir -p ~/workshop/gsea
cd ~/workshop/gsea
mkdir -p fgsea_example genabel real_data_gsea results
```

Your structure should now look like:

```text
~/workshop/
├── gwas_*.assoc.logistic
├── srd_qc_allchr.*
├── srd_qc_hwe.*
└── gsea/
    ├── fgsea_example/
    ├── genabel/
    ├── real_data_gsea/
    └── results/
```

---

<details>
<summary><strong>Shared GSEA resources</strong></summary>

<br>

For the complete real-data exercise, the workshop server contains the following shared resources:

```text
/workshop/data/gsea/
├── all_pathways.gmt
├── srd_snp_gene_map.tsv
└── ARS-UCD1.2_gene.csv
```

Only the first two are required for the main student workflow.

### `srd_snp_gene_map.tsv`

This should contain at least:

```text
SNP    GENE
```

Each row links one SNP from the workshop dataset to a gene.

A SNP may appear more than once if it maps to more than one gene.

### `all_pathways.gmt`

This contains the biological gene sets used for GSEA in standard GMT format.

The pathway collection used in this tutorial includes gene sets from:

- BioCarta
- Gene Ontology
- KEGG
- PANTHER
- Reactome

If the pathway collection is updated, the gene identifiers in the pathway file must match the identifiers in the SNP-to-gene mapping file.

### Optional annotation file

`ARS-UCD1.2_gene.csv` can be retained if students will also learn how the SNP-to-gene mapping file was generated.

</details>

---

# Part 1 — In-class GSEA using `fgsea`

## 1. What is GSEA?

**Gene Set Enrichment Analysis (GSEA)** asks whether members of a predefined group of genes tend to occur disproportionately near the top of a ranked gene list.

Instead of asking only:

> Which individual SNPs reached genome-wide significance?

GSEA asks:

> Are multiple genes belonging to the same biological process collectively showing stronger association signals than expected?

This makes GSEA useful when a biological pathway contains several genes with moderate effects that may not individually pass a stringent GWAA significance threshold.

### GSEA versus GSEA-SNP

Traditional GSEA was developed for gene-expression data.

**GSEA-SNP** extends the idea to genome-wide association data:

```text
GWAA SNPs
   ↓
SNP-level association statistics
   ↓
SNPs assigned to genes
   ↓
Genes ranked by association evidence
   ↓
Gene sets/pathways tested for enrichment
```

---

## 2. Start the `fgsea` example

Move into the example directory:

```bash
cd ~/workshop/gsea/fgsea_example
```

Start R:

```bash
R
```

Load the packages:

```r
library(fgsea)
library(data.table)
library(ggplot2)
```

If `fgsea` does not load, click below.

<details>
<summary><strong>Package not found? Click here.</strong></summary>

<br>

`fgsea` is installed through **Bioconductor**.

If package installation is permitted on your account:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("fgsea")
```

Then load it:

```r
library(fgsea)
```

> ⚠️ On the workshop server, packages should normally already be installed. Do not repeatedly reinstall packages if `library(fgsea)` works.

</details>

---

## 3. Load the built-in example data

The `fgsea` package contains example ranked genes and pathways.

Run:

```r
data(examplePathways)
data(exampleRanks)

set.seed(12345)
```

Inspect them:

```r
head(exampleRanks)
length(exampleRanks)

head(examplePathways)
length(examplePathways)
```

### What are we looking at?

`exampleRanks` is a **named numeric vector**:

```text
Gene ID  -> ranking statistic
```

`examplePathways` is a **list of pathways**:

```text
Pathway name -> genes belonging to that pathway
```

Try:

```r
head(names(exampleRanks))
head(names(examplePathways))
```

Then inspect one pathway:

```r
examplePathways[[1]]
```

### Question 🤔

What is the difference between:

```r
names(exampleRanks)
```

and:

```r
examplePathways[[1]]
```

<details>
<summary><strong>Click to reveal the answer</strong></summary>

<br>

`names(exampleRanks)` gives the identifiers of the genes in the **ranked gene list**.

`examplePathways[[1]]` gives the genes belonging to **one predefined biological gene set**.

GSEA determines whether genes belonging to a pathway occur unusually high in the ranked list.

</details>

---

## 4. Run the example GSEA

To match the in-class example, use:

```r
fgsea_results <- fgsea(
    pathways = examplePathways,
    stats    = exampleRanks,
    minSize  = 1,
    maxSize  = 5000
)
```

Inspect the results:

```r
head(fgsea_results)
```

Sort them by adjusted p-value:

```r
fgsea_results <- fgsea_results[order(padj)]

head(fgsea_results, 10)
```

### Main result columns

| Column | Meaning |
|---|---|
| `pathway` | Name of the pathway/gene set |
| `pval` | Enrichment p-value |
| `padj` | Multiple-testing-adjusted p-value |
| `ES` | Enrichment score |
| `NES` | Normalized enrichment score |
| `size` | Number of pathway genes represented in the ranked list |
| `leadingEdge` | Genes contributing most strongly to the enrichment signal |

> 💡 The broad `minSize = 1` and `maxSize = 5000` settings are used here for the classroom example. In a research analysis, biologically sensible pathway-size filters such as 10–500 or 15–500 are often preferable.

---

## 5. What is the enrichment score?

Imagine walking down the ranked list from the strongest gene to the weakest gene.

- When a gene belongs to the pathway, the running score increases.
- When a gene does not belong to the pathway, the score decreases.
- The largest deviation from zero is the **enrichment score (ES)**.

The **normalized enrichment score (NES)** adjusts the enrichment score so that pathways of different sizes can be compared more fairly.

### Leading-edge genes

The **leading-edge subset** contains the genes contributing most strongly to the pathway enrichment.

These are especially useful when interpreting a significant pathway.

---

## 6. Plot one pathway

The in-class tutorial uses the programmed cell death pathway:

```r
p_one <- plotEnrichment(
    examplePathways[["5991130_Programmed_Cell_Death"]],
    exampleRanks
) +
    labs(title = "Programmed Cell Death")

p_one
```

Save it:

```r
ggsave(
    "fgsea_programmed_cell_death.png",
    plot = p_one,
    width = 8,
    height = 5,
    dpi = 300
)
```

---

## 7. Plot several pathways

Identify strongly enriched pathways:

```r
topPathwaysUp <- fgsea_results[
    ES > 0
][head(order(pval), n = 10), pathway]

topPathwaysDown <- fgsea_results[
    ES < 0
][head(order(pval), n = 10), pathway]

topPathways <- c(
    topPathwaysUp,
    rev(topPathwaysDown)
)
```

Plot them:

```r
plotGseaTable(
    examplePathways[topPathways],
    exampleRanks,
    fgsea_results,
    gseaParam = 0.5
)
```

Save the table plot:

```r
png(
    "fgsea_top_pathways.png",
    width = 4800,
    height = 1800,
    res = 300
)

plotGseaTable(
    examplePathways[topPathways],
    exampleRanks,
    fgsea_results,
    gseaParam = 0.5
)

dev.off()
```

---

## 8. Export the pathway results needed for homework

We will save **two types of results**:

- **Full results:** all pathways tested by `fgsea`, useful for the checkpoint and comparisons.
- **Filtered homework results:** only pathways with **`NES > 2.5` and `padj < 0.05`**, including all pathway genes represented in the ranked list and the leading-edge genes.

> **Important:** The homework filter is for **strong positively enriched pathways**. It does **not** include negatively enriched pathways, even if their `padj` is below 0.05. Use the **unfiltered `fgsea_results`** when answering general questions about all significant pathways or the largest absolute NES.

### 8A. Save the complete result table (CSV)

The `leadingEdge` column is a list, so convert it to text in a **copy**, not in the original `fgsea_results` object.

```r
fgsea_export <- copy(fgsea_results)

fgsea_export[, leadingEdge := vapply(
    leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
)]

fwrite(fgsea_export, "fgsea_example_results.csv")
```

### 8B. Prepare the filtered homework table

Copy and run this block **while you are still in the R session** and while `fgsea_results`, `examplePathways`, and `exampleRanks` are available.

```r
# 1. Filter significant pathways (strong positive enrichment)
fgsea_export_sig <- copy(
  fgsea_results[NES > 2.5 & padj < 0.05]
)

# 2. Sort by decreasing NES
setorder(fgsea_export_sig, -NES)

# 3. Extract all pathway genes used in GSEA
fgsea_export_sig[, all_genes := lapply(
  pathway,
  function(p) {
    intersect(
      examplePathways[[p]],
      names(exampleRanks)
    )
  }
)]

# 4. Count leading-edge genes
fgsea_export_sig[, leading_edge_count := lengths(leadingEdge)]

# 5. Verify pathway sizes
stopifnot(
  all(lengths(fgsea_export_sig$all_genes) ==
        fgsea_export_sig$size)
)

# 6. Convert gene lists to semicolon-separated strings
fgsea_export_sig[, all_genes :=
  vapply(all_genes, paste, collapse = ";", FUN.VALUE = character(1))
]

fgsea_export_sig[, leadingEdge :=
  vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
]

# 7. Export filtered CSV
csv_file <- "fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.csv"
fwrite(fgsea_export_sig, csv_file)

# 8. Export the same table as Excel (.xlsx)
xlsx_file <- "fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.xlsx"

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(as.data.frame(fgsea_export_sig), xlsx_file)
  message("Excel file saved: ", xlsx_file)
} else {
  message("The CSV was saved, but writexl is not installed; ",
          "see the Excel note below.")
}

# Confirm how many pathways were exported
cat("Filtered pathways exported:", nrow(fgsea_export_sig), "\n")
```

The exported table includes the original GSEA statistics plus:

| Column | Meaning |
|---|---|
| `all_genes` | All genes in that pathway that were present in `exampleRanks` and used in GSEA |
| `leadingEdge` | The leading-edge subset from the enrichment analysis |
| `leading_edge_count` | Number of genes in that leading-edge subset |

> **Why the semicolons?** CSV and Excel cells cannot directly store an R list column in a useful way. Separating gene IDs with `;` keeps the full set of IDs in one cell. **Do not replace the original list columns in `fgsea_results`**, since they may be needed for later plots and comparisons.

<details>
<summary><strong>Excel file not generated? Optional package setup</strong></summary>

If the workbook is missing, you can **use the CSV file for homework** or open it in Excel and save it as `.xlsx` with **File → Save As**. CSV is the primary file and contains the same data.

If your R environment allows package installation, install the writer once:

```r
install.packages("writexl")
```

Then repeat only the Excel export:

```r
writexl::write_xlsx(
  as.data.frame(fgsea_export_sig),
  "fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.xlsx"
)
```

On the workshop server, **do not spend class time troubleshooting package installation**; the CSV output is sufficient.

</details>

### 8C. Verify the exports

```r
list.files(pattern = "^fgsea_example.*\\.(csv|xlsx)$")
```

When finished, exit R:

```r
q("no")
```

In the terminal (in `~/workshop/gsea/fgsea_example`), confirm:

```bash
ls -lh fgsea_example*.csv fgsea_example*.xlsx 2>/dev/null
```

The expected files are:

```text
fgsea_example_results.csv
fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.csv
fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.xlsx  (if writexl is installed)
fgsea_programmed_cell_death.png
fgsea_top_pathways.png
```

**For homework, use the filtered CSV or Excel workbook** to examine the pathways, their genes, and leading-edge genes. Keep the full CSV as well: it lets you answer questions involving all pathways, including those that do not pass the homework filter.

---

## 🧠 Checkpoint

Answer these before continuing:

1. How many pathways have `padj < 0.05`?
2. Which pathway has the largest absolute NES?
3. What genes make up the leading-edge subset of the top pathway?
4. Why do we use `padj` rather than looking only at `pval`?
5. What information is contained in the ranked list that is not contained in the pathway file?


> 🛑 **STOP HERE FOR CLASS.** Once you have exported the homework results and answered the five questions, **do not proceed to GenABEL or the real-data GSEA workflow during class**. Those sections are optional and available for future exploration in [02_Optional_Project.md](02_Optional_Project.md).

**Homework files to keep:** `fgsea_example_sig_results_nes_2.5_all_genes_leading_edge.csv` and, if generated, the matching `.xlsx` file.


---

[Return to the README](README.md) · [Optional project](02_Optional_Project.md)
