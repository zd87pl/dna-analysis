# Canonical GRCh38 coordinates used by Helixight

Single source of truth for the variant positions hard-coded in `scripts/*.sh`,
documented in `docs/variant_database.md` and shipped in `data/sample_genome.vcf`.
All positions are **GRCh38**, 1-based, on **bare contig names** (`1`..`22`, `X`,
`Y`, `MT`; `scripts/lib/vcf_common.sh` renames `chr`-prefixed input before any
query). `REF>ALT` is given as it appears in a plus-strand VCF record; for genes
on the minus strand this is the complement of the cDNA/HGVS allele (for example
MTHFR C677T is `G>A` on the plus strand).

Confidence column:

- `canonical` - row from the synthesizer's canonical table; scripts, docs and the
  sample fixture must agree with it.
- `unverified` - recollected without a dbSNP lookup (the review container has no
  network access). Keep, but leave a `# TODO(build): unverified GRCh38` comment
  next to any script literal that relies on it, and check against dbSNP before
  relying on it clinically.

The `hx_gt_allele REGION REF ALT [RSID]` check in the shared library refuses to
score a record whose REF/ALT (or rsID) does not match the row, so a wrong
coordinate degrades to "not reported" plus a stderr warning rather than a wrong
call.

## Canonical table

| Gene | rsID | GRCh38 chrom:pos | REF>ALT | Interpretation note | Confidence |
|------|------|------------------|---------|---------------------|------------|
| ACTN3 | rs1815739 | 11:66560624 | C>T | R577X; T/T = X/X (no alpha-actinin-3), C/C = R/R power-associated | canonical |
| ACE | rs4646994 (alias rs1799752) | 17:63488529 | anchor+287 bp Alu > anchor | Alu insertion/deletion: REF carries the insertion (I), ALT is the deletion (D); 0/0 = I/I, 0/1 = I/D, 1/1 = D/D. Query with `hx_gt_indel 17:63488529 250` | canonical |
| COMT | rs4680 | 22:19963748 | G>A | Val158Met; A = Met (worrier, lower enzyme activity). Was 22:19951271 (GRCh37) | canonical |
| CYP1A2 | rs762551 | 15:74749576 | A>C | A = *1F fast caffeine metaboliser; C = slow (non-*1F). Was 15:75041917 | canonical |
| APOE | rs429358 | 19:44908684 | T>C | C defines the e4 allele (with rs7412 C). Was 19:45411941 | canonical |
| APOE | rs7412 | 19:44908822 | C>T | T defines the e2 allele (with rs429358 T). Was 19:45412079 | canonical |
| LCT/MCM6 | rs4988235 | 2:135851076 | G>A | -13910C>T; A (T on cDNA) = lactase persistence. Was 2:136608646 | canonical |
| ALDH2 | rs671 | 12:111803962 | G>A | ALDH2*2; A = alcohol flush, reduced acetaldehyde clearance. Was 12:112241766 | canonical |
| FTO | rs9939609 | 16:53786615 | T>A | A = obesity risk allele. Was 16:53820527 | canonical |
| TCF7L2 | rs7903146 | 10:112998590 | C>T | T = type 2 diabetes risk allele. Was 10:114758349 | canonical |
| HFE | rs1800562 | 6:26092913 | G>A | C282Y; A/A = highest hemochromatosis risk. Was 6:26093141 (and 6:26091179 in carrier_screening, which is H63D on GRCh37) | canonical |
| HFE | rs1799945 | 6:26090951 | C>G | H63D; mild risk modifier, compound heterozygote with C282Y matters | canonical |
| MTHFR | rs1801133 | 1:11796321 | G>A | C677T (minus strand); A = 677T, reduced enzyme activity. Was 1:11856378 | canonical |
| MTHFR | rs1801131 | 1:11794419 | T>G | A1298C (minus strand); G = 1298C. Was 1:11854476 | canonical |
| AGT | rs699 | 1:230710048 | A>G | M235T; G = 235T, higher angiotensinogen. Was 1:230845794 | canonical |
| BDNF | rs6265 | 11:27658369 | C>T | Val66Met (minus strand); T = Met. Was 11:27679916 | canonical |
| CRP | rs1205 | 1:159712443 | C>T | 3'UTR (3872G>A in gene orientation, minus strand); T (minor) = lower baseline CRP, C/C = highest. Was 1:159682233 | canonical |
| SOD2 | rs4880 | 6:159692840 | A>G | Val16Ala (minus strand); G = Ala. Was 6:160113872 | canonical |
| GDF5 | rs143383 | 20:35437977 | A>G (allele orientation unverified) | 5'UTR variant, osteoarthritis / height association; A = OA risk allele (+104T on the minus-strand cDNA), G = protective. Was 20:34025756 | canonical (position) |
| HERC2 | rs12913832 | 15:28120472 | A>G | G/G = blue eyes most likely, A/A = brown. Was 15:28365618 | canonical |
| VDR | rs1544410 | 12:47846052 | C>T | BsmI (minus strand) | canonical |
| VDR | rs2228570 | 12:47879112 | A>G | FokI (minus strand); scripts used GRCh37 12:48272895 mislabelled as rs1544410 | canonical |
| PPARG | rs1801282 | 3:12351626 | C>G | Pro12Ala; G = Ala, improved insulin sensitivity. PRS used 6:20679709 | canonical |
| CYP2D6 | rs3892097 | 22:42128945 | C>T | *4 (splice defect); T = no function. pharmacogenomics used 22:42126611 (= rs1135840) | canonical |
| CYP2D6 | rs1065852 | 22:42130692 | G>A | *10 (P34S); A = decreased function | canonical |
| CYP2D6 | rs16947 | 22:42127941 | G>A | *2 (R296C); normal function marker | canonical |
| SLCO1B1 | rs4149056 | 12:21178615 | T>C | *5 (V174A); C = reduced transport, simvastatin myopathy risk. Was 12:21176804 (= rs2306283) | canonical |
| TPMT | rs1800462 | 6:18143724 | C>G | *2 (A80P); G = reduced function | canonical |
| TPMT | rs1800460 | 6:18138997 | C>T | *3B (A154T); T = reduced function | canonical |
| TPMT | rs1142345 | 6:18130687 | T>C | *3C (Y240C); C = reduced function | canonical |
| GBA | rs76763715 | 1:155235843 | T>C | N370S (N409S); C = Gaucher allele | canonical |
| GBA | rs421016 | 1:155235252 | A>G | L444P (L483P); G = Gaucher allele | canonical |
| ABCC11 | rs17822931 | 16:48224287 | C>T | C/C and C/T = wet earwax, T/T = dry. Was 16:48258198 | canonical |
| ADRB2 | rs1042713 | 5:148826877 | G>A | Arg16Gly; already GRCh38 in scripts | canonical |
| PPARGC1A | rs8192678 | 4:23814039 | C>T | Gly482Ser (minus strand); T = Ser; already GRCh38 in scripts | canonical |
| IL6 | rs1800795 | 7:22727026 | G>C | -174G/C promoter; already GRCh38 in scripts | canonical |
| CDKN2B-AS1 | rs1333049 | 9:22125504 | G>C | 9p21 coronary artery disease risk locus | unverified |
| CDKN2B-AS1 | rs10811661 | 9:22134095 | T>C | 9p21 type 2 diabetes risk locus | unverified |
| SLC30A8 | rs13266634 | 8:117172544 | C>T | R325W; C = type 2 diabetes risk allele | unverified |
| ZEB2 region | rs10427255 | ~2:145.37 Mb (exact position unknown) | C>T | Photic sneeze reflex. NOT 2:135851076 (that is rs4988235). Dropped from the sample fixture until verified | unverified |

## Additional rows carried by the sample fixture

Rows below are not part of the canonical table; they are positions the scripts
already query on GRCh38 or well-known ClinVar/PharmGKB variants. They are
listed so that `data/sample_genome.vcf` has one documented source. Treat them
as `unverified` for clinical purposes.

| Gene | rsID | GRCh38 chrom:pos | REF>ALT | Interpretation note | Confidence |
|------|------|------------------|---------|---------------------|------------|
| ABO | rs8176719 | 9:133257521 (or left-aligned 133257520) | T>TC | 261delG: the reference carries the O (deletion) allele, so the ALT insertion = non-O; 0/0 = O/O, 0/1 = X/O, 1/1 = non-O/non-O. blood_type.sh also accepts the deletion form TC>T (ALT = O). Query with `hx_gt_indel 9:133257520-133257521 1` and match POS/REF/ALT, not with a SNP lookup | unverified (orientation/POS) |
| ABO | rs8176746 | 9:133256264 | G>T | B allele marker (blood_type.sh); c.796C>A on the minus strand, so the plus-strand VCF form is G>T: G = non-B (A/O), T = B allele | unverified |
| ABO | rs7853989 | 9:133256042 | G>C | A vs B marker (blood_type.sh) | unverified |
| RHD (proxy) | rs590787 | 1:25362249 | G>A | Rh proxy position queried by blood_type.sh; rsID/position not independently verified | unverified |
| RHD (proxy) | rs676785 | 1:25409008 | C>T | Rh proxy position queried by blood_type.sh; rsID/position not independently verified | unverified |
| CFTR | rs113993960 | 7:117559590 | ATCT>A | F508del, 3-bp in-frame deletion (ClinVar VCF representation) | unverified |
| HBB | rs334 | 11:5227002 | T>A | HbS sickle allele | unverified |
| HBB | rs33930165 | 11:5227003 | C>T | HbC (c.19G>A, the base next to HbS on the minus strand); derived from the rs334 anchor (carrier_screening.sh) | unverified |
| HBB | rs33971440 | 11:5226929 | C>T | Beta-thalassemia IVS-I-1 c.92+1G>A (the G>T form is C>A); derived from the rs334 anchor via exon 1 = 92 coding nt (carrier_screening.sh) | unverified |
| HEXA | rs76173977 | 15:72346579 | T>TTATC | Tay-Sachs c.1274_1277dupTATC, 4-bp insertion; carrier_screening.sh and mega_analysis.sh (the former mega_analysis literal 15:72638892 was the GRCh37 locus) | unverified |
| PPARA | rs4253778 | 22:46150500 | G>C | G = endurance-associated allele (Ahmetov 2006); one literal shared by athletic_genetics.sh and triathlon_genetics.sh, queried with hx_gt_allele. The value is probably the GRCh38 gene start rather than the intronic SNP (~22:46219976 expected) and the former athletic literal 22:46615880 was GRCh37-range; verify against dbSNP before relying on it | unverified |
| PAH | rs5030858 | 12:102840632 | G>A | R408W | unverified |
| F5 | rs6025 | 1:169549811 | C>T | Factor V Leiden (minus strand) | unverified |
| CYP2C19 | rs4244285 | 10:94781859 | G>A | *2 loss of function | unverified |
| CYP2C19 | rs4986893 | 10:94780653 | G>A | *3 loss of function | unverified |
| CYP2C19 | rs12248560 | 10:94761900 | C>T | *17 increased activity | unverified |
| CYP2C9 | rs1799853 | 10:94942290 | C>T | *2 | unverified |
| CYP2C9 | rs1057910 | 10:94981296 | A>C | *3 | unverified |
| CYP3A4 | rs35599367 | 7:99768693 | C>T | *22 reduced expression | unverified |
| CYP3A5 | rs776746 | 7:99672916 | T>C | *3 splice defect (C = non-expresser) | unverified |
| VKORC1 | rs9923231 | 16:31096368 | C>T | -1639G>A (minus strand); T = warfarin sensitivity. pharmacogenomics.sh lists it as G>A (cDNA orientation) | unverified |
| DPYD | rs3918290 | 1:97915614 | C>T | *2A splice variant | unverified |
| OPRM1 | rs1799971 | 6:154039662 | A>G | A118G | unverified |
| ADH1B | rs1229984 | 4:99318162 | T>C | His48Arg: REF T = His48 (*2, fast acetaldehyde production); ALT C = Arg48 (*1, common European allele); queried by mega_analysis.sh and personalized_protocol.sh | unverified |
| IRF4 | rs12203592 | 6:396321 | C>T | Freckling / sun sensitivity; the old sample carried this position under the MC1R rsID rs1805007 | unverified |
| MC1R | rs1805007 | 16:89919436 | C>T | R151C red-hair allele; matches docs/variant_database.md and mega_analysis.sh (old sample position 16:89736156 is outside MC1R) | unverified |
| TAS2R38 | rs713598 | 7:141972804 | C>G | A49P PTC taster marker; position queried by mega_analysis.sh (old sample labelled it rs72921001/OR6A2) | unverified |
| MTNR1B | rs10830963 | 11:92975544 | C>G | Melatonin receptor / fasting glucose; old sample had 15:74219546 | unverified |
