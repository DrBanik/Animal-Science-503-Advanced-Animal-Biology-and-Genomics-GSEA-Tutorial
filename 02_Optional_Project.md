# 02 — Optional Project: GenABEL and GSEA on Bovine GWAA Results

> **OPTIONAL — NOT REQUIRED FOR THE IN-CLASS SESSION.** This is a reference for students who want to explore the full workflow after class. **The required in-class material ends at the checkpoint in [01_In_Class_Tutorial.md](01_In_Class_Tutorial.md).**

This file preserves the advanced parts of the original Neibergs Lab–adapted tutorial: GenABEL permutations; additive GWAA to SNP-to-gene ranking; `fgsea` analysis on bovine data; interpretation, QC, and troubleshooting. **Each major section is collapsed by default**. Click a section title to expand it whenever you want to explore.

Before starting, complete the `fgsea` example from [01_In_Class_Tutorial.md](01_In_Class_Tutorial.md), and confirm that the required GWAA files and shared gene/pathway resources are available on the server. For this optional workflow, use **additive GWAA results with valid SNP association statistics**; not every covariate combination will necessarily produce usable results.

[Return to the README](README.md) · [In-class tutorial](01_In_Class_Tutorial.md)

---

<details>
<summary><strong>Part 2 — GenABEL (click to expand — optional)</strong></summary>

# Part 2 — GenABEL

## Why are we using GenABEL?

In this tutorial, **GenABEL** is used to generate phenotype-permuted genome-wide association statistics.

The major steps were:

```text
PED + MAP + phenotype
          ↓
       GenABEL
          ↓
       IBS matrix
          ↓
    Polygenic model
          ↓
Environmental residuals
          ↓
Phenotype permutations
          ↓
Permuted SNP statistics
```

These permutation statistics provide an empirical null distribution of SNP-level association signals.

The workshop uses `fgsea` for the final enrichment analysis, so the GenABEL permutation files are **not direct inputs to `fgsea`**. This section is included to show how phenotype permutations are generated and why an empirical null distribution can be useful in GSEA-SNP analyses.

> ⚠️ `GenABEL` is an older package. The workshop server should provide a compatible R environment. Do not spend workshop time trying to install GenABEL into a modern R environment unless instructed to do so.

---

## 9. Prepare the GWAA dataset for GenABEL

Move into the GenABEL directory:

```bash
cd ~/workshop/gsea/genabel
```

The GWAA tutorial produced a binary PLINK dataset:

```text
~/workshop/srd_qc_allchr.bed
~/workshop/srd_qc_allchr.bim
~/workshop/srd_qc_allchr.fam
```

GenABEL can read PLINK-style PED/MAP data, so convert the binary files:

```bash
plink \
    --bfile ~/workshop/srd_qc_allchr \
    --autosome-num 30 \
    --allow-no-sex \
    --recode \
    --out srd_gsea
```

### Main outputs

```text
srd_gsea.ped
srd_gsea.map
```

Check them:

```bash
head srd_gsea.map
head -n 2 srd_gsea.ped
wc -l srd_gsea.map
wc -l srd_gsea.ped
```

---

## 10. Create the GenABEL phenotype file

GenABEL requires a phenotype file containing an `id` column and a `sex` column.

For this workshop, we will create it directly from the filtered PLINK `.fam` file.

Create an R script:

```bash
nano make_genabel_pheno.R
```

Paste:

```r
fam <- read.table(
    "../../srd_qc_allchr.fam",
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(fam) <- c(
    "FID",
    "IID",
    "Father",
    "Mother",
    "Sex_PLINK",
    "Phenotype_PLINK"
)

genabel_pheno <- data.frame(
    id = fam$IID,

    # GenABEL sex coding convention used here:
    # female = 0, male = 1
    sex = ifelse(
        fam$Sex_PLINK == 1, 1,
        ifelse(fam$Sex_PLINK == 2, 0, 0)
    ),

    # PLINK binary phenotype:
    # 1 = control, 2 = case
    # Convert to:
    # 0 = control, 1 = case
    status = ifelse(
        fam$Phenotype_PLINK == 2, 1,
        ifelse(fam$Phenotype_PLINK == 1, 0, NA)
    )
)

if (any(is.na(genabel_pheno$status))) {
    stop("Missing or unexpected phenotype values were found.")
}

write.table(
    genabel_pheno,
    "genabel_pheno.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

print(head(genabel_pheno))
print(table(genabel_pheno$status))
```

Save with **Ctrl+O**, press **Enter**, then exit with **Ctrl+X**.

Run:

```bash
Rscript make_genabel_pheno.R
```

Inspect the file:

```bash
head genabel_pheno.txt
```

You should see:

```text
id    sex    status
...
```

---

## 11. Start GenABEL

Start R, then load GenABEL:

```r
library(GenABEL)
```

Confirm that the package loaded:

```r
packageVersion("GenABEL")
```

Now navigate to the GenABEL directory if needed:

```r
setwd("~/workshop/gsea/genabel")
```

---

## 12. Save the marker names

Create a marker-name file so the permutation statistics can always be matched back to the correct SNP.

```r
marker_map <- read.table(
    "srd_gsea.map",
    header = FALSE,
    stringsAsFactors = FALSE
)

marker.names <- marker_map$V2

write.table(
    marker.names,
    "marker.names.txt",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
```

Check:

```r
head(marker.names)
length(marker.names)
```

The number of marker names should match the number of SNPs in your MAP file.

---

## 13. Convert PED/MAP data to GenABEL format

Run:

```r
convert.snp.ped(
    pedfile = "srd_gsea.ped",
    mapfile = "srd_gsea.map",
    outfile = "srd_gsea.raw",
    strand = "u",
    mapHasHeaderLine = FALSE
)
```

### What does this do?

- `pedfile` specifies the genotype/phenotype PED file.
- `mapfile` specifies the SNP map.
- `outfile` is the GenABEL-formatted genotype file.
- `strand = "u"` means the strand is treated as unknown.
- `mapHasHeaderLine = FALSE` tells GenABEL that the PLINK MAP file has no header.

---

## 14. Load genotype and phenotype data together

```r
srd_data <- load.gwaa.data(
    phe = "genabel_pheno.txt",
    gen = "srd_gsea.raw",
    force = TRUE
)
```

Inspect the phenotype data:

```r
head(phdata(srd_data))
```

Check the number of SNPs and animals:

```r
nsnps(srd_data)
nids(srd_data)
```

### Question 🤔

Do these counts agree with the final filtered GWAA dataset?

If they do not, stop and determine where the mismatch occurred before continuing.

---

## 15. Create the IBS matrix

Account for genomic similarity among animals using an **identity-by-state (IBS)** matrix.

```r
srd_gkin <- ibs(
    srd_data[, autosomal(srd_data)],
    weight = "freq"
)
```

Save it:

```r
saveRDS(
    srd_gkin,
    "srd_gkin.rds"
)

write.table(
    srd_gkin,
    "srd_gkin.txt",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
)
```

> 💡 Saving the matrix prevents you from having to recalculate it if a later step needs to be rerun.

---

## 16. Fit the polygenic model

For this teaching example, use the binary `status` variable:

```r
poly_srd <- polygenic(
    status,
    kinship.matrix = srd_gkin,
    data = srd_data
)
```

Use the environmental residuals from the polygenic model:

```r
newtrait <- as.numeric(
    poly_srd$pgresidualY
)
```

Check:

```r
summary(newtrait)
length(newtrait)
```

### Why use `pgresidualY`?

These residuals represent the portion of the phenotype remaining after the polygenic component has been accounted for.

These residuals are used to generate phenotype permutations after the polygenic component has been accounted for. Shuffling the residuals breaks the phenotype-genotype correspondence while leaving the genotype data and relatedness structure unchanged.

---

## 17. Demonstrate phenotype permutations

A full permutation analysis can use a much larger number of permutations:

```text
1000 permutations per file
10 output files
```

That is too computationally expensive for a classroom demonstration.

For the workshop, we will generate only **10 permutations**.

Add a temporary phenotype column:

```r
srd_data@phdata$perm_trait <- NA_real_
```

Run:

```r
n.perms <- 10

perm_chi2 <- matrix(
    NA_real_,
    nrow = length(marker.names),
    ncol = n.perms
)

for (i in seq_len(n.perms)) {

    cat(
        "Permutation",
        i,
        "of",
        n.perms,
        "\n"
    )

    srd_data@phdata$perm_trait <- sample(
        newtrait,
        replace = FALSE
    )

    qtscore_perm <- qtscore(
        perm_trait,
        data = srd_data,
        clambda = FALSE
    )

    p_perm <- as.vector(
        qtscore_perm[, "Pc1df"]
    )

    perm_chi2[, i] <- qchisq(
        p_perm,
        df = 1,
        lower.tail = FALSE
    )
}
```

Attach marker names:

```r
perm_demo <- data.frame(
    Marker = marker.names,
    perm_chi2,
    check.names = FALSE
)

colnames(perm_demo)[-1] <- paste0(
    "CHI2_PERM_",
    seq_len(n.perms)
)
```

Save:

```r
write.table(
    perm_demo,
    "genabel_permutation_demo.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

Exit R when finished:

```r
q("no")
```

Inspect the result:

```bash
head genabel_permutation_demo.tsv
```

### Main GenABEL outputs

```text
genabel_pheno.txt
marker.names.txt
srd_gsea.raw
srd_gkin.rds
srd_gkin.txt
genabel_permutation_demo.tsv
```

---

<details>
<summary><strong>What does a full GenABEL permutation analysis look like?</strong></summary>

<br>

A full analysis can generate much larger null distributions, for example:

```r
n.perms <- 1000
n.files <- 10
```

Each permutation randomized the polygenic environmental residual and reran the SNP score test.

The resulting chi-square statistics represented SNP association signals expected when the phenotype-genotype relationship had been broken by permutation.

The resulting empirical null distribution can then be used to compare the observed pathway enrichment against what is expected after permutation.

The classroom demonstration uses only 10 permutations because the purpose here is to understand the logic rather than reproduce a computationally expensive production run.

</details>

---

</details>

---

<details>
<summary><strong>Part 3 — Run GSEA on your GWAA results (click to expand — optional)</strong></summary>

# Part 3 — Run GSEA on your GWAA results

Now we will leave the example data behind and use the **same GWAA results you generated in the previous tutorial**.

For consistency, we will use an **additive GWAA model**.

---

## 18. Choose the additive GWAA result

Move into the real-data directory:

```bash
cd ~/workshop/gsea/real_data_gsea
```

The GWAA tutorial generated several additive models.

List them:

```bash
ls ~/workshop/gwas_*_ADD.assoc.logistic
```

For a reproducible classroom example, we will use:

```text
gwas_unadjusted_ADD.assoc.logistic
```

Create a symbolic link so we do not duplicate the file:

```bash
ln -sf \
    ~/workshop/gwas_unadjusted_ADD.assoc.logistic \
    gwas_additive.assoc.logistic
```

Check:

```bash
ls -lh gwas_additive.assoc.logistic
```

> 💡 If you want to use a different covariate-adjusted additive model, change only the source file above. For example, `gwas_sire_ADD.assoc.logistic` can be linked instead. The downstream GSEA workflow remains the same.

---

## 19. Inspect the additive GWAA result

Run:

```bash
head gwas_additive.assoc.logistic
```

The PLINK logistic association output should contain columns such as:

```text
CHR
SNP
BP
A1
TEST
NMISS
OR
STAT
P
```

We only want rows corresponding to the additive SNP test:

```text
TEST = ADD
```

### Important

GSEA uses **all valid SNPs**, not only genome-wide significant SNPs.

Do **not** filter the GWAA results to `p < 0.05`, `p < 1e-5`, or genome-wide significance before creating the ranked list.

---

## 20. Convert additive SNP p-values to chi-square statistics

Create:

```bash
nano prepare_gsea_ranks.R
```

Paste:

```r
library(data.table)

assoc <- fread(
    "gwas_additive.assoc.logistic"
)

assoc <- assoc[
    TEST == "ADD" &
    !is.na(P) &
    P > 0 &
    P <= 1
]

assoc[, CHI2 := qchisq(
    P,
    df = 1,
    lower.tail = FALSE
)]

observed_snp_stats <- assoc[
    ,
    .(
        SNP,
        CHR,
        BP,
        P,
        CHI2
    )
]

fwrite(
    observed_snp_stats,
    "observed_additive_snp_chi2.tsv",
    sep = "\t"
)

cat(
    "SNPs retained:",
    nrow(observed_snp_stats),
    "\n"
)

print(
    observed_snp_stats[
        order(-CHI2)
    ][1:10]
)
```

Run:

```bash
Rscript prepare_gsea_ranks.R
```

Inspect:

```bash
head observed_additive_snp_chi2.tsv
```

### Why chi-square?

For this tutorial, additive GWAA p-values are converted to 1-degree-of-freedom chi-square statistics.

This produces an **unsigned association-strength statistic**:

```text
larger chi-square = stronger SNP-level evidence
```

---

## 21. Map SNPs to genes

For the workshop, use the shared SNP-to-gene mapping file:

```text
/workshop/data/gsea/srd_snp_gene_map.tsv
```

Inspect it:

```bash
head /workshop/data/gsea/srd_snp_gene_map.tsv
```

The required columns are:

```text
SNP    GENE
```

### Important

The SNP-to-gene mapping should include **all mapped SNPs**, not only significant SNPs.

This is essential because GSEA evaluates the distribution of association evidence across the complete ranked list.

---

<details>
<summary><strong>How was the SNP-to-gene mapping created?</strong></summary>

<br>

The SNP-to-gene mapping can be created by expanding gene coordinates by an analysis-specific genomic window and assigning SNPs that fall inside those regions to genes.

A modern implementation can perform the same overlap operation using genomic-range software such as `GenomicRanges`.

The exact window should be defined by the study design and should **not** be chosen after looking at which value gives the most significant pathways.

For this workshop, the mapping is supplied so everyone uses the same SNP-to-gene definition.

</details>

---

## 22. Convert SNP statistics to gene-level statistics

A gene can contain or be assigned multiple SNPs.

To create one value per gene, we will use the **largest chi-square statistic among SNPs assigned to that gene**.

This mirrors the idea that the SNP with the strongest evidence represents the gene in a GSEA-SNP ranking.

Create:

```bash
nano make_gene_ranks.R
```

Paste:

```r
library(data.table)

snp_stats <- fread(
    "observed_additive_snp_chi2.tsv"
)

snp_gene <- fread(
    "/workshop/data/gsea/srd_snp_gene_map.tsv"
)

required_cols <- c(
    "SNP",
    "GENE"
)

if (!all(required_cols %in% names(snp_gene))) {
    stop(
        "SNP-to-gene file must contain columns named SNP and GENE."
    )
}

mapped <- merge(
    snp_gene,
    snp_stats,
    by = "SNP",
    all = FALSE
)

mapped <- mapped[
    !is.na(GENE) &
    GENE != ""
]

gene_stats <- mapped[
    ,
    {
        best <- which.max(CHI2)

        .(
            STAT = CHI2[best],
            TOP_SNP = SNP[best],
            TOP_SNP_P = P[best],
            N_MAPPED_SNPS = .N
        )
    },
    by = GENE
]

gene_stats <- gene_stats[
    order(-STAT)
]

fwrite(
    gene_stats,
    "gene_rank_statistics.tsv",
    sep = "\t"
)

cat(
    "Mapped SNPs:",
    nrow(mapped),
    "\n"
)

cat(
    "Genes in ranked list:",
    nrow(gene_stats),
    "\n"
)

print(
    head(gene_stats, 10)
)
```

Run:

```bash
Rscript make_gene_ranks.R
```

Inspect:

```bash
head gene_rank_statistics.tsv
```

---

## 🧠 Checkpoint

Before running GSEA:

```bash
wc -l observed_additive_snp_chi2.tsv
wc -l /workshop/data/gsea/srd_snp_gene_map.tsv
wc -l gene_rank_statistics.tsv
```

### Questions

1. Why are there fewer genes than SNPs?
2. Why can one gene have several mapped SNPs?
3. Why are we retaining only the strongest SNP statistic for each gene?
4. What could happen if genes with many SNPs systematically have more chances to obtain an extreme statistic?

<details>
<summary><strong>Think about Question 4 before clicking 😏</strong></summary>

<br>

Genes represented by many SNPs have more opportunities to contain an extreme association statistic.

This is a potential **gene-size/SNP-density bias** when the maximum SNP statistic is used to represent a gene.

Permutation-based GSEA-SNP approaches can help account for properties of the observed SNP structure through empirical permutation. In this workshop, the simplified `fgsea` analysis is intended to teach the pathway-enrichment workflow and should be interpreted with this limitation in mind.

</details>

---

## 23. Inspect the pathway file

The workshop pathway file is:

```text
/workshop/data/gsea/all_pathways.gmt
```

Inspect it:

```bash
head -n 3 /workshop/data/gsea/all_pathways.gmt
```

A GMT file generally contains:

```text
Pathway_Name    Description    Gene1    Gene2    Gene3 ...
```

The identifiers must match the `GENE` identifiers in:

```text
gene_rank_statistics.tsv
```

---

## 24. Run `fgsea` on the GWAA-derived gene ranking

Create:

```bash
nano run_real_data_fgsea.R
```

Paste:

```r
library(data.table)
library(fgsea)
library(ggplot2)

set.seed(12345)

# ------------------------------------------------------------
# 1. Read gene-level ranking
# ------------------------------------------------------------

gene_stats <- fread(
    "gene_rank_statistics.tsv"
)

ranks <- gene_stats$STAT
names(ranks) <- gene_stats$GENE

ranks <- sort(
    ranks,
    decreasing = TRUE
)

# Remove any non-finite values
ranks <- ranks[
    is.finite(ranks)
]

# ------------------------------------------------------------
# 2. Read pathway definitions
# ------------------------------------------------------------

pathways <- gmtPathways(
    "/workshop/data/gsea/all_pathways.gmt"
)

cat(
    "Genes in ranked list:",
    length(ranks),
    "\n"
)

cat(
    "Pathways loaded:",
    length(pathways),
    "\n"
)

# ------------------------------------------------------------
# 3. Run GSEA
# ------------------------------------------------------------

fgsea_results <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    scoreType = "pos"
)

fgsea_results <- fgsea_results[
    order(padj, -NES)
]

# ------------------------------------------------------------
# 4. Export results
# ------------------------------------------------------------

fgsea_export <- copy(
    fgsea_results
)

fgsea_export[, leadingEdge := vapply(
    leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
)]

fwrite(
    fgsea_export,
    "../results/SRD_additive_GSEA_results.csv"
)

# ------------------------------------------------------------
# 5. Export significant pathways
# ------------------------------------------------------------

sig_results <- fgsea_export[
    !is.na(padj) &
    padj < 0.05
]

fwrite(
    sig_results,
    "../results/SRD_additive_GSEA_FDR05.csv"
)

# ------------------------------------------------------------
# 6. Print top results
# ------------------------------------------------------------

cat(
    "\nNumber of pathways with FDR < 0.05:",
    nrow(sig_results),
    "\n\n"
)

print(
    head(
        fgsea_export[
            ,
            .(
                pathway,
                pval,
                padj,
                ES,
                NES,
                size
            )
        ],
        15
    )
)

# ------------------------------------------------------------
# 7. Plot the top pathway
# ------------------------------------------------------------

if (nrow(fgsea_results) > 0) {

    top_pathway <- fgsea_results$pathway[1]

    p_top <- plotEnrichment(
        pathways[[top_pathway]],
        ranks
    ) +
        labs(
            title = top_pathway
        )

    ggsave(
        "../results/SRD_top_pathway_enrichment.png",
        plot = p_top,
        width = 9,
        height = 5.5,
        dpi = 300
    )
}

# ------------------------------------------------------------
# 8. Plot a table of the top pathways
# ------------------------------------------------------------

n_show <- min(
    10,
    nrow(fgsea_results)
)

if (n_show > 0) {

    top_pathways <- fgsea_results$pathway[
        seq_len(n_show)
    ]

    png(
        "../results/SRD_top_pathways_table.png",
        width = 2600,
        height = 2000,
        res = 300
    )

    plotGseaTable(
        pathways[top_pathways],
        ranks,
        fgsea_results,
        gseaParam = 0.5
    )

    dev.off()
}
```

Run:

```bash
Rscript run_real_data_fgsea.R
```

---

## 25. Why did we use `scoreType = "pos"`?

Look back at the ranking statistic:

```text
CHI2
```

A chi-square statistic cannot be negative.

Our ranking therefore measures:

```text
weak association  -------------------->  strong association
```

rather than:

```text
protective direction  <---->  risk direction
```

Because the statistics are unsigned and we are interested only in enrichment toward the **strong-association end** of the ranked list, we use:

```r
scoreType = "pos"
```

This differs from the built-in `fgsea` example, where the ranking statistics can contain both positive and negative values.

---

## 26. Inspect the results

Return to Bash and inspect:

```bash
cd ~/workshop/gsea/results
```

List:

```bash
ls -lh
```

You should have:

```text
SRD_additive_GSEA_results.csv
SRD_additive_GSEA_FDR05.csv
SRD_top_pathway_enrichment.png
SRD_top_pathways_table.png
```

View the first few results:

```bash
head SRD_additive_GSEA_results.csv
```

If X11 is working, view the figures:

```bash
feh SRD_top_pathway_enrichment.png
```

and:

```bash
feh SRD_top_pathways_table.png
```

If graphical forwarding is unavailable, download the images using `scp` or FileZilla as you did in the GWAA tutorial.

---

</details>

---

<details>
<summary><strong>Part 4 — Interpreting your GSEA results (click to expand — optional)</strong></summary>

# Part 4 — Interpreting your GSEA results

## 27. Adjusted p-value

The main significance column is:

```text
padj
```

For this exercise, use:

```text
FDR < 0.05
```

as the primary multiple-testing threshold.

Count significant pathways in Bash:

```bash
awk -F',' 'NR > 1 && $3 < 0.05 {n++} END {print n}' \
    SRD_additive_GSEA_results.csv
```

Or inspect them directly in R:

```r
results <- data.table::fread(
    "SRD_additive_GSEA_results.csv"
)

results[
    padj < 0.05
]
```

---

## 28. NES

The **normalized enrichment score (NES)** reflects the strength of pathway enrichment after accounting for gene-set size and the GSEA null distribution.

For this analysis, larger positive NES values indicate that genes belonging to the pathway are concentrated toward the **stronger GWAA-association end** of the ranked gene list.

Do not interpret the NES as an effect size for fetal loss.

It is a **pathway enrichment statistic**.

---

## 29. Leading-edge genes

The most biologically informative genes are often contained in:

```text
leadingEdge
```

These genes drive much of the enrichment signal.

Inspect the top pathway:

```r
results <- data.table::fread(
    "SRD_additive_GSEA_results.csv"
)

results[
    order(padj)
][1, .(
    pathway,
    padj,
    NES,
    leadingEdge
)]
```

### Next question

Do any of the leading-edge genes contain or lie near SNPs that were among the strongest signals in your GWAA?

That comparison connects the pathway-level result back to the Manhattan plot.

---

## 30. Trace a leading-edge gene back to its SNP

Suppose a leading-edge gene is called:

```text
GENE_X
```

Find it in the gene-ranking file:

```bash
grep -w "GENE_X" \
    ~/workshop/gsea/real_data_gsea/gene_rank_statistics.tsv
```

The output tells you:

```text
GENE
STAT
TOP_SNP
TOP_SNP_P
N_MAPPED_SNPS
```

Now find the SNP in the additive GWAA:

```bash
grep -w "TOP_SNP_NAME" \
    ~/workshop/gwas_unadjusted_ADD.assoc.logistic
```

This lets you move from:

```text
Pathway
   ↓
Leading-edge gene
   ↓
Top SNP assigned to gene
   ↓
GWAA association
```

That connection is important when biologically interpreting GSEA results.

---

</details>

---

<details>
<summary><strong>Part 5 — Questions for the workshop (click to expand — optional)</strong></summary>

# Part 5 — Questions for the workshop

## GSEA concepts

1. Why can a pathway be significantly enriched even if none of its individual SNPs reaches genome-wide significance?
2. What is the difference between ES and NES?
3. Why is `padj` more appropriate than raw `pval` when many pathways are tested?
4. What is a leading-edge gene?
5. Why should pathway size be restricted with `minSize` and `maxSize`?

## GWAA → GSEA

6. Why did we use all valid additive-model SNPs instead of only significant SNPs?
7. Why did we convert p-values to chi-square statistics?
8. Why did we use `scoreType = "pos"` for the real GWAA analysis?
9. What information is lost when an unsigned chi-square statistic is used?
10. Why might genes containing many SNPs be advantaged when the maximum SNP statistic is used?

## Interpretation

11. Which pathway has the smallest FDR-adjusted p-value?
12. Which pathway has the largest NES?
13. How many pathways have FDR < 0.05?
14. Which genes appear in the leading-edge subset of the top pathway?
15. Are any of those genes also located near strong signals in the Manhattan plot?
16. Do multiple significant pathways share the same leading-edge genes?
17. If several pathways contain many of the same genes, should they automatically be interpreted as independent biological findings?

---

</details>

---

<details>
<summary><strong>Part 6 — Useful quality-control checks (click to expand — optional)</strong></summary>

# Part 6 — Useful quality-control checks

## Check 1: Number of additive SNP results

```bash
awk '$5 == "ADD" {n++} END {print n}' \
    ~/workshop/gwas_unadjusted_ADD.assoc.logistic
```

## Check 2: Number of SNPs with valid chi-square statistics

```bash
wc -l \
    ~/workshop/gsea/real_data_gsea/observed_additive_snp_chi2.tsv
```

Remember that the file contains one header row.

## Check 3: Number of genes

```bash
wc -l \
    ~/workshop/gsea/real_data_gsea/gene_rank_statistics.tsv
```

## Check 4: Number of pathways tested

In R:

```r
length(pathways)
nrow(fgsea_results)
```

These numbers can differ because pathways smaller than `minSize`, larger than `maxSize`, or with too little overlap with the ranked gene list are excluded.

## Check 5: Pathway-to-ranked-list overlap

```r
length(
    intersect(
        pathways[[1]],
        names(ranks)
    )
)
```

---

</details>

---

<details>
<summary><strong>Part 7 — Common problems (click to expand — optional)</strong></summary>

# Part 7 — Common problems

<details>
<summary><strong>Problem: fgsea says all statistics are positive</strong></summary>

<br>

That is expected for this workshop because the ranking statistic is chi-square.

Make sure the real-data analysis includes:

```r
scoreType = "pos"
```

</details>

---

<details>
<summary><strong>Problem: Very few pathways are tested</strong></summary>

<br>

Check whether the identifiers in the pathway file match the identifiers in the ranked gene list.

In R:

```r
head(names(ranks))
head(pathways[[1]])
```

Then calculate overlap:

```r
length(
    intersect(
        names(ranks),
        unique(unlist(pathways))
    )
)
```

If overlap is extremely low, the pathway file and SNP-to-gene mapping may use different gene identifier systems.

</details>

---

<details>
<summary><strong>Problem: My gene-ranking file is almost empty</strong></summary>

<br>

Check the SNP identifiers:

```bash
head observed_additive_snp_chi2.tsv
head /workshop/data/gsea/srd_snp_gene_map.tsv
```

Then compare:

```r
library(data.table)

a <- fread(
    "observed_additive_snp_chi2.tsv"
)

m <- fread(
    "/workshop/data/gsea/srd_snp_gene_map.tsv"
)

length(
    intersect(
        a$SNP,
        m$SNP
    )
)
```

If the overlap is very small, the SNP-to-gene mapping file does not correspond to the marker naming system used by the GWAA dataset.

</details>

---

<details>
<summary><strong>Problem: GenABEL will not load</strong></summary>

<br>

Do **not** install random versions of GenABEL during the workshop.

GenABEL is an older package, but it is already configured on the workshop server.

Confirm:

```r
R.version.string
```

and:

```r
.libPaths()
```

If `library(GenABEL)` still fails, stop here and check the workshop R setup before continuing.

</details>

---

</details>

---

<details>
<summary><strong>Part 8 — Quick command summary (click to expand — optional)</strong></summary>

# Part 8 — Quick command summary

| Step | Main command/script | Main output |
|---|---|---|
| Login | `ssh -Y studentXX@10.104.58.24` | Server session |
| Example GSEA | `fgsea(...)` | See the CSV/XLSX exports in `01_In_Class_Tutorial.md` |
| Prepare GenABEL data | `plink --recode` | `srd_gsea.ped`, `srd_gsea.map` |
| GenABEL phenotype | `Rscript make_genabel_pheno.R` | `genabel_pheno.txt` |
| GenABEL IBS | `ibs(...)` | `srd_gkin.*` |
| GenABEL permutations | `qtscore(...)` loop | `genabel_permutation_demo.tsv` |
| Additive GWAA input | `gwas_*_ADD.assoc.logistic` | SNP p-values |
| Convert p → χ² | `Rscript prepare_gsea_ranks.R` | `observed_additive_snp_chi2.tsv` |
| SNP → gene | `Rscript make_gene_ranks.R` | `gene_rank_statistics.tsv` |
| Real-data GSEA | `Rscript run_real_data_fgsea.R` | GSEA result CSVs + plots |
| Interpret | inspect `padj`, `NES`, `leadingEdge` | Candidate pathways/genes |

---


---

</details>

---

# References and further reading

1. Subramanian A, et al. (2005). **Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.** *Proceedings of the National Academy of Sciences* 102(43):15545–15550.  
   DOI: 10.1073/pnas.0506580102

2. Holden M, Deng S, Wojnowski L, Kulle B. (2008). **GSEA-SNP: applying gene set enrichment analysis to SNP data from genome-wide association studies.** *Bioinformatics* 24:2784–2785.  
   DOI: 10.1093/bioinformatics/btn516

3. Aulchenko YS, Ripke S, Isaacs A, van Duijn CM. (2007). **GenABEL: an R library for genome-wide association analysis.** *Bioinformatics* 23(10):1294–1296.

4. `fgsea` Bioconductor tutorial:  
   https://bioconductor.org/packages/release/bioc/html/fgsea.html

---

---

# Notes for students ✍️📖

- Keep the **GWAA and GSEA directories separate**, but do not delete the GWAA outputs.
- GSEA depends on the SNP identifiers remaining consistent across the GWAA results, SNP-to-gene mapping, and genotype data.
- Always inspect files with `head`, `wc -l`, and `ls -lh` before running a long analysis.
- Do not filter the GWAA results to significant SNPs before GSEA.
- Read the first R error carefully before rerunning a script.
- A significant pathway does **not** mean that every gene in that pathway is associated with the phenotype.
- Pay special attention to the **leading-edge genes**.
- Multiple pathway databases contain overlapping biological information, so related significant pathways are not necessarily independent findings.
- Treat pathway enrichment as a tool for biological interpretation and hypothesis generation alongside the GWAA—not as a replacement for the GWAA.

---

## The End! 🧬🐄🎉

At this point you should be able to move from:

```text
GWAA SNP association results
```

to:

```text
biological pathways and leading-edge genes
```

while understanding how phenotype permutation with GenABEL relates to pathway enrichment with `fgsea`.
