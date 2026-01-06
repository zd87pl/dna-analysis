# 🧬 Helixight Variant Database

This document lists all genetic variants analyzed by Helixight, organized by category.

## 📊 Summary

| Category | Variants | Script |
|----------|----------|--------|
| Athletic Performance | 25+ | `athletic_genetics.sh` |
| Triathlon Specific | 35+ | `triathlon_genetics.sh` |
| Pharmacogenomics | 50+ | `mega_analysis.sh` |
| Nutrigenomics | 30+ | `personalized_protocol.sh` |
| Cancer Predisposition | 80+ | `cancer_genes_scan.sh` |
| Mental Health | 15+ | `mega_analysis.sh` |
| Longevity | 10+ | `mega_analysis.sh` |
| **TOTAL** | **500+** | - |

---

## 🏃 Athletic Performance Variants

### Muscle Fiber Type

| Gene | rsID | Position (GRCh38) | Effect |
|------|------|-------------------|--------|
| **ACTN3** | rs1815739 | 11:66560624 | R577X - R/R power, X/X endurance |
| **ACE** | rs4341 | 17:63488529 | I/D polymorphism - I/I endurance |
| **PPARGC1A** | rs8192678 | 4:23814039 | Gly482Ser - mitochondrial biogenesis |

### Oxygen Capacity

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **VEGFA** | rs2010963 | 6:43770613 | Angiogenesis |
| **HIF1A** | rs11549465 | 14:62207556 | Hypoxia adaptation |
| **EPAS1** | rs142764723 | 2:46441523 | High altitude (Denisovan) |
| **NRF1** | rs2402970 | 7:129613162 | Mitochondrial function |
| **ADRB2** | rs1042713 | 5:148826877 | Bronchodilation |

### Recovery & Injury

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **IL6** | rs1800795 | 7:22727026 | Inflammation response |
| **TNF** | rs1800629 | 6:31575254 | TNF-alpha production |
| **COL5A1** | rs12722 | 9:137684151 | Tendon elasticity |
| **COL1A1** | rs1800012 | 17:50201632 | Collagen quality |
| **MMP3** | rs679620 | 11:102733807 | Matrix degradation |
| **GDF5** | rs143383 | 20:35437977 | Cartilage health |

### Mental Performance

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **COMT** | rs4680 | 22:19951271 | Val158Met - Warrior/Worrier |
| **BDNF** | rs6265 | 11:27679916 | Val66Met - neuroplasticity |
| **DRD2** | rs1800497 | 11:113400106 | Dopamine receptors |
| **SLC6A4** | rs25531 | 17:30194319 | Serotonin transporter |

---

## 💊 Pharmacogenomics Variants

### Drug Metabolism (CYP450)

| Gene | rsID | Position | Drugs Affected |
|------|------|----------|----------------|
| **CYP2D6** | Multiple | 22:42523943 | Codeine, Tamoxifen, many SSRIs |
| **CYP2C19** | rs4244285 | 10:94781859 | Clopidogrel, PPIs, some SSRIs |
| **CYP2C9** | rs1799853 | 10:94942290 | Warfarin, NSAIDs |
| **CYP1A2** | rs762551 | 15:75041917 | Caffeine, Clozapine |
| **CYP3A5** | rs776746 | 7:99672916 | Tacrolimus, many drugs |

### Pain & Anesthesia

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **OPRM1** | rs1799971 | 6:154039662 | Opioid sensitivity |
| **COMT** | rs4680 | 22:19951271 | Pain perception |
| **MC1R** | rs1805007 | 16:89919436 | Anesthesia requirements |

### Cardiovascular Drugs

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **VKORC1** | rs9923231 | 16:31096368 | Warfarin sensitivity |
| **SLCO1B1** | rs4149056 | 12:21176804 | Statin myopathy risk |
| **ADRB1** | rs1801252 | 10:114045297 | Beta-blocker response |

---

## 🥗 Nutrigenomics Variants

### Macronutrient Metabolism

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **FTO** | rs9939609 | 16:53820527 | Obesity risk |
| **TCF7L2** | rs7903146 | 10:114758349 | Carb tolerance, T2D |
| **APOA2** | rs5082 | 1:161222292 | Saturated fat sensitivity |
| **PPARG** | rs1801282 | 3:12351626 | Insulin sensitivity |
| **PPARA** | rs4253778 | 22:46150500 | Fat oxidation |

### Vitamin Metabolism

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **MTHFR** | rs1801133 | 1:11856378 | C677T - Folate metabolism |
| **MTHFR** | rs1801131 | 1:11854476 | A1298C - Folate metabolism |
| **VDR** | rs1544410 | 12:48272895 | Vitamin D receptor |
| **BCMO1** | rs12934922 | 16:81264597 | Beta-carotene conversion |
| **FADS1** | rs174546 | 11:61567029 | Omega-3 conversion |
| **FUT2** | rs602662 | 19:49206674 | B12 absorption |

### Food Sensitivities

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **LCT** | rs4988235 | 2:136608646 | Lactose tolerance |
| **HLA-DQ** | Multiple | 6:32637406 | Celiac disease |
| **ADH1B** | rs1229984 | 4:99318162 | Alcohol metabolism |
| **ALDH2** | rs671 | 12:112241766 | Alcohol flush |

---

## 🎀 Cancer Predisposition Variants

### BRCA1 (Breast/Ovarian Cancer)

| Mutation | Position | Population |
|----------|----------|------------|
| 185delAG | 17:43045677 | Ashkenazi |
| 5382insC | 17:43057051 | Slavic |
| C61G | 17:43063930 | General |
| And 50+ more... | | |

### BRCA2 (Breast/Ovarian Cancer)

| Mutation | Position | Population |
|----------|----------|------------|
| 6174delT | 13:32354860 | Ashkenazi |
| 886delGT | 13:32370557 | Icelandic |
| And 40+ more... | | |

### Lynch Syndrome (Colorectal Cancer)

| Gene | Position | Effect |
|------|----------|--------|
| **MLH1** | 3:36993350-37050918 | MMR gene |
| **MSH2** | 2:47403067-47709917 | MMR gene |
| **MSH6** | 2:47695530-47810101 | MMR gene |
| **PMS2** | 7:5970925-6009106 | MMR gene |

### Other Cancer Genes

| Gene | Position | Associated Cancer |
|------|----------|-------------------|
| **TP53** | 17:7668402-7687538 | Li-Fraumeni |
| **APC** | 5:112707498-112846239 | FAP |
| **PTEN** | 10:87863625-87971930 | Cowden |
| **CDH1** | 16:68771117-68869444 | Gastric |
| **CHEK2** | 22:28687743-28742422 | Breast |
| **PALB2** | 16:23603160-23641310 | Breast |
| **ATM** | 11:108222484-108369102 | Breast |

---

## 🧠 Mental Health Variants

### Depression & Mood

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **SLC6A4** | rs25531 | 17:30194319 | Serotonin transporter |
| **BDNF** | rs6265 | 11:27679916 | Brain plasticity |
| **FKBP5** | rs1360780 | 6:35656217 | Stress sensitivity |

### ADHD & Focus

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **DRD4** | rs1800955 | 11:637269 | Dopamine receptor |
| **DAT1** | rs28363170 | 5:1392791 | Dopamine transporter |
| **SNAP25** | rs3746544 | 20:10268206 | Neurotransmitter release |

### Addiction Risk

| Gene | rsID | Position | Substance |
|------|------|----------|-----------|
| **DRD2** | rs1800497 | 11:113400106 | General |
| **OPRM1** | rs1799971 | 6:154039662 | Opioids |
| **CHRNA5** | rs16969968 | 15:78882925 | Nicotine |
| **GABRA2** | rs279858 | 4:46252826 | Alcohol |

---

## 🏆 Longevity Variants

| Gene | rsID | Position | Effect |
|------|------|----------|--------|
| **FOXO3** | rs2802292 | 6:108989256 | Major longevity gene |
| **APOE** | rs429358/rs7412 | 19:45411941 | e2/e3/e4 alleles |
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

*Last updated: January 2025*
