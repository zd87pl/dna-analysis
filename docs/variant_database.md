# 🧬 Helixight Variant Database

This document lists the genetic variants analyzed by Helixight, organized by category.

All positions are **GRCh38**, 1-based, on bare contig names (`1`..`22`, `X`, `Y`, `MT`).
The rows that the scripts score are kept in sync with the canonical table in
[`scripts/lib/GRCh38_COORDINATES.md`](../scripts/lib/GRCh38_COORDINATES.md); positions
that are not in that table were taken from the scripts and have not been re-verified
against dbSNP, so treat them as informational.

## 📊 Summary

The "Markers" column is the number of markers each script queries, i.e. the denominator
of the `HX_DATA_COVERAGE: found/queried markers` line every script prints at the end
(measured on `data/sample_genome.vcf.gz`). Many markers are shared between scripts.

| Category | Markers | Script |
|----------|---------|--------|
| Athletic Performance | 18 | `athletic_genetics.sh` |
| Triathlon Specific | 30 | `triathlon_genetics.sh` |
| Personalized Protocol (training, supplements, nutrigenomics) | 26 | `personalized_protocol.sh` |
| Mega Analysis (10 sections: ancestry, carrier status, mental health, chronotype, nutrigenomics, immunity, senses, skin/hair/ageing, longevity, supplementation) | 91 checks at 81 positions | `mega_analysis.sh` |
| Fun Genetics (traits, ancestry hints) | 35 | `fun_genetics.sh` |
| Pharmacogenomics | 26 | `pharmacogenomics.sh` |
| Carrier Screening | 29 | `carrier_screening.sh` |
| Blood Type Prediction | 5 | `blood_type.sh` |
| Polygenic Risk Scores | 31 | `polygenic_risk_scores.sh` |
| Cancer Predisposition | 13 gene regions + 17 point probes (14 named BRCA1/BRCA2/TP53/Lynch mutations, CHEK2 1100delC at 22:28695868, and two position-only MUTYH probes at 1:45332879 and 1:45340175) | `cancer_genes_scan.sh` |
| Rare Variants | 24 clinical gene regions (no fixed markers) | `rare_variants.sh` |
| ClinVar Scan | every Pathogenic/Likely_pathogenic record of the downloaded ClinVar release | `clinvar_scan.sh` |
| mtDNA Haplogroup | MT contig extract for Haplogrep (no fixed markers) | `extract_mtdna_haplogrep.sh` |
| **TOTAL** | **291 single-marker checks** (18+30+26+91+35+26+29+5+31) **+ 17 mutation sites + 37 gene regions** | - |

Across all scripts the repository holds fewer than 400 distinct variant keys (rsIDs and
chrom:pos literals combined), so counts above that are not supported by the code.

---

## 🏃 Athletic Performance Variants

### Muscle Fiber Type

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **ACTN3** | rs1815739 | 11:66560624 | R577X - R/R power, X/X endurance |
| **ACE** | rs4646994 (alias rs1799752) | 17:63488529 | Alu I/D polymorphism - I/I endurance. In the VCF the REF allele carries the 287 bp Alu insertion (I) and ALT is the deletion (D): 0/0 = I/I, 0/1 = I/D, 1/1 = D/D |
| **PPARGC1A** | rs8192678 | 4:23814039 | Gly482Ser - mitochondrial biogenesis |

### Oxygen Capacity

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **VEGFA** | rs2010963 | 6:43770613 | Angiogenesis |
| **HIF1A** | rs11549465 | 14:62207556 | Hypoxia adaptation |
| **EPAS1** | rs142764723 | 2:46441523 | High altitude (Denisovan) |
| **NRF1** | rs2402970 | 7:129613162 | Mitochondrial function |
| **ADRB2** | rs1042713 | 5:148826877 | Bronchodilation |

### Recovery & Injury

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **IL6** | rs1800795 | 7:22727026 | Inflammation response |
| **TNF** | rs1800629 | 6:31575254 | TNF-alpha production |
| **COL5A1** | rs12722 | not probed (former 9:137684151 is outside GRCh38 COL5A1, 9:134.64-134.84 Mb; unverified) | Tendon elasticity |
| **COL1A1** | rs1800012 | 17:50201632 | Collagen quality |
| **MMP3** | rs679620 | 11:102733807 | Matrix degradation |
| **GDF5** | rs143383 | 20:35437977 | Cartilage health |

### Mental Performance

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **COMT** | rs4680 | 22:19963748 | Val158Met - Warrior/Worrier |
| **BDNF** | rs6265 | 11:27658369 | Val66Met - neuroplasticity |
| **DRD2** | rs1800497 | 11:113400106 | Dopamine receptors |
| **SLC6A4** | rs25531 | 17:30194319 | Serotonin transporter |

---

## 💊 Pharmacogenomics Variants

### Drug Metabolism (CYP450)

| Gene | rsID | Position (GRCh38) | Drugs Affected |
|------|------|-------------------|----------------|
| **CYP2D6** | rs3892097 (*4), rs1065852 (*10), rs16947 (*2) | 22:42128945, 22:42130692, 22:42127941 | Codeine, Tamoxifen, many SSRIs |
| **CYP2C19** | rs4244285 | 10:94781859 | Clopidogrel, PPIs, some SSRIs |
| **CYP2C9** | rs1799853 | 10:94942290 | Warfarin, NSAIDs |
| **CYP1A2** | rs762551 | 15:74749576 | Caffeine, Clozapine |
| **CYP3A5** | rs776746 | 7:99672916 | Tacrolimus, many drugs |

### Pain & Anesthesia

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **OPRM1** | rs1799971 | 6:154039662 | Opioid sensitivity |
| **COMT** | rs4680 | 22:19963748 | Pain perception |
| **MC1R** | rs1805007 | 16:89919436 | Anesthesia requirements |

### Cardiovascular Drugs

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **VKORC1** | rs9923231 | 16:31096368 | Warfarin sensitivity |
| **SLCO1B1** | rs4149056 | 12:21178615 | Statin myopathy risk |
| **ADRB1** | rs1801252 | 10:114045297 | Beta-blocker response |

---

## 🥗 Nutrigenomics Variants

### Macronutrient Metabolism

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **FTO** | rs9939609 | 16:53786615 | Obesity risk |
| **TCF7L2** | rs7903146 | 10:112998590 | Carb tolerance, T2D |
| **APOA2** | rs5082 | 1:161222292 | Saturated fat sensitivity |
| **PPARG** | rs1801282 | 3:12351626 | Insulin sensitivity |
| **PPARA** | rs4253778 | 22:46150500 (unverified; G = endurance-associated allele) | Fat oxidation |

### Vitamin Metabolism

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **MTHFR** | rs1801133 | 1:11796321 | C677T - Folate metabolism |
| **MTHFR** | rs1801131 | 1:11794419 | A1298C - Folate metabolism |
| **VDR** | rs1544410 | 12:47846052 | Vitamin D receptor (BsmI) |
| **VDR** | rs2228570 | 12:47879112 | Vitamin D receptor (FokI) |
| **BCMO1** | rs12934922 | 16:81264597 | Beta-carotene conversion |
| **FADS1** | rs174546 | 11:61567029 | Omega-3 conversion |
| **FUT2** | rs602662 | 19:49206674 | B12 absorption |

### Food Sensitivities

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **LCT/MCM6** | rs4988235 | 2:135851076 | Lactose tolerance |
| **HLA-DQ** | Multiple | 6:32637406 | Celiac disease |
| **ADH1B** | rs1229984 | 4:99318162 | Alcohol metabolism |
| **ALDH2** | rs671 | 12:111803962 | Alcohol flush |

---

## 🎀 Cancer Predisposition Variants

### BRCA1 (Breast/Ovarian Cancer)

The BRCA1/BRCA2 point-mutation positions below are **unverified** GRCh38
coordinates derived from the standard exon layout and the GRCh37->GRCh38
offsets implied by the scan's gene regions (BRCA1 +1,847,983; BRCA2
-574,137). Because a position-only match cannot confirm the named mutation,
the scan reports an indel-shaped record at a `del`/`ins` position (a SNP at
M1775R) only as an informational note; these hits are **not** counted as
pathogenic findings. Use the ClinVar scan for a verified call.

| Mutation | Position (GRCh38, derived) | Population |
|----------|----------------------------|------------|
| 185delAG (c.68_69delAG) | ~17:43124027 | Ashkenazi |
| 5382insC (c.5266dupC) | ~17:43057062 | Slavic |
| M1775R (c.5324T>G) | ~17:43051071 | BRCT domain |

### BRCA2 (Breast/Ovarian Cancer)

| Mutation | Position (GRCh38, derived) | Population |
|----------|----------------------------|------------|
| 6174delT (c.5946delT) | ~13:32340300 | Ashkenazi |
| 4075delGT (c.3847_3848delGT) | ~13:32338202 | Pathogenic |
| 886delGT (c.658_659delGT) | ~13:32329469 | Icelandic |

### CHEK2 / MUTYH

| Gene | Variant | Position (GRCh38) | Match |
|------|---------|-------------------|-------|
| CHEK2 | 1100delC | 22:28695868 (unverified) | indel records only |
| MUTYH | recessive, counted together | 1:45332879, 1:45340175 (unverified) | SNP records only, position-only |

### TP53 (Li-Fraumeni) hotspots

| Mutation | Position (GRCh38) | Population |
|----------|-------------------|------------|
| R175H | 17:7675088 (C>T, rs28934578) | Hotspot |
| R248Q | 17:7674220 (C>T, rs11540652) | Hotspot |
| R248W | 17:7674221 (G>A, rs121912651) | Hotspot |
| R273H | 17:7673802 (C>T, rs28934576) | Hotspot |
| R282W | 17:7673776 (G>A, rs28934574) | Pathogenic |

Positions are derived from the GRCh37 dbSNP coordinates (+96,682, the offset of
the TP53 gene region) and are matched by REF/ALT/rsID; they are marked
*unverified* in `scripts/cancer_genes_scan.sh` until checked against ClinVar.
(An earlier table listed 17:7676154 as "R282W"; that position is rs1042522
Pro72Arg, a benign polymorphism, and is no longer probed.)

### Lynch Syndrome (Colorectal Cancer)

| Gene | Position (GRCh38) | Effect |
|------|-------------------|--------|
| **MLH1** | 3:36993350-37050918 | MMR gene (listed mutation: c.350C>T at 3:37042337) |
| **MSH2** | 2:47403067-47709917 | MMR gene (listed mutation: c.1216C>T at 2:47641560) |
| **MSH6** | 2:47695530-47810101 | MMR gene (listed mutation: c.3261dup at 2:47703303) |
| **PMS2** | 7:5970925-6009106 | MMR gene |

### Other Cancer Genes

| Gene | Position (GRCh38) | Associated Cancer |
|------|-------------------|-------------------|
| **TP53** | 17:7668402-7687538 | Li-Fraumeni |
| **APC** | 5:112707498-112846239 | FAP |
| **PTEN** | 10:87863625-87971930 | Cowden |
| **CDH1** | 16:68737292-68835537 | Gastric |
| **CHEK2** | 22:28687743-28742422 | Breast |
| **PALB2** | 16:23603160-23641310 | Breast |
| **ATM** | 11:108222484-108369102 | Breast |

---

## 🧠 Mental Health Variants

### Depression & Mood

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **SLC6A4** | rs25531 | 17:30194319 | Serotonin transporter |
| **BDNF** | rs6265 | 11:27658369 | Brain plasticity |
| **FKBP5** | rs1360780 | 6:35656217 | Stress sensitivity |

### ADHD & Focus

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **DRD4** | rs1800955 | 11:637269 | Dopamine receptor |
| **DAT1** | rs28363170 | 5:1392791 | Dopamine transporter |
| **SNAP25** | rs3746544 | 20:10268206 | Neurotransmitter release |

### Addiction Risk

| Gene | rsID | Position (GRCh38) | Substance |
|------|------|-------------------|-----------|
| **DRD2** | rs1800497 | 11:113400106 | General |
| **OPRM1** | rs1799971 | 6:154039662 | Opioids |
| **CHRNA5** | rs16969968 | 15:78882925 | Nicotine |
| **GABRA2** | rs279858 | 4:46252826 | Alcohol |

---

## 🏆 Longevity Variants

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **FOXO3** | rs2802292 | 6:108989256 | Major longevity gene |
| **APOE** | rs429358 / rs7412 | 19:44908684 / 19:44908822 | e2/e3/e4 alleles |
| **CETP** | rs5882 | 16:57015091 | HDL cholesterol |
| **TERT** | rs2736100 | 5:1287194 | Telomere length |
| **IL6** | rs1800795 | 7:22727026 | Inflammation |

---

## 📚 References

1. NCBI ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
2. PharmGKB: https://www.pharmgkb.org/
3. SNPedia: https://www.snpedia.com/
4. gnomAD: https://gnomad.broadinstitute.org/
5. GWAS Catalog: https://www.ebi.ac.uk/gwas/

---

*Last updated: September 2026 (Helixight 1.3.0)*
