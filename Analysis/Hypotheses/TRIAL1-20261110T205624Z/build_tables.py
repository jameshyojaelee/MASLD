"""Builds trial_outcomes.tsv and excluded.tsv for TRIAL1. Every value below was copied from a
PubMed abstract (NCBI efetch) or a ClinicalTrials.gov v2 record opened on 2026-11-10."""
import csv, os
D = os.path.dirname(os.path.abspath(__file__))
NV = "not verified"
cols = ["mechanism","drug","target_genes","action_direction","chembl_id","nct_id","phase","n_randomized",
        "duration_weeks","population","endpoint","result_vs_control","effect_size","ci_or_p","pmid","doi",
        "verification","notes","representative"]
R = []
def row(*a): R.append(dict(zip(cols, a)))

ABS, REG, AR = "PubMed abstract", "ClinicalTrials.gov registry record (results section)", "PubMed abstract; ClinicalTrials.gov registry record"

# ---- THR-beta
m=("THR-beta agonist","resmetirom","THRB","activate","CHEMBL3261331")
row(*m,"NCT03900429","3","966 (primary analysis population; registry enrollment 1759)","52","biopsy NASH, fibrosis F1B-F3","mash_resolution","improved","100 mg 29.9% vs placebo 9.7% (80 mg 25.9%)","P<0.001 both doses","38324483","10.1056/NEJMoa2309000",AR,"MAESTRO-NASH; dual primary endpoints","yes")
row(*m,"NCT03900429","3","966 (primary analysis population; registry enrollment 1759)","52","biopsy NASH, fibrosis F1B-F3","fibrosis_improvement","improved","100 mg 25.9% vs placebo 14.2% (80 mg 24.2%)","P<0.001 both doses","38324483","10.1056/NEJMoa2309000",AR,"MAESTRO-NASH; dual primary endpoints","yes")
row(*m,"NCT04197479","3","972 double-blind (325/327/320) plus 171 open-label","52","NAFLD with presumed NASH (noninvasive)","pdff","improved","16-week hepatic fat LSM difference vs placebo: 100 mg -38.6%, 80 mg -34.9%; 52-week -33.9, -28.8","P<0.0001 (16 week)","37845512","10.1038/s41591-023-02603-1",AR,"MAESTRO-NAFLD-1; key secondary endpoint; abstract says 'hepatic fat' and does not name the MRI method","yes")
# ---- GLP-1R
m=("GLP-1R agonist","semaglutide","GLP1R","activate","CHEMBL2108724")
row(*m,"NCT04822181","3","1197 (interim analysis of first 800)","72","biopsy MASH, fibrosis F2-F3","mash_resolution","improved","62.9% vs 34.3%; difference 28.7 percentage points","95% CI 21.1 to 36.2; P<0.001","40305708","10.1056/NEJMoa2413258",AR,"ESSENCE part 1; 2.4 mg weekly","yes")
row(*m,"NCT04822181","3","1197 (interim analysis of first 800)","72","biopsy MASH, fibrosis F2-F3","fibrosis_improvement","improved","36.8% vs 22.4%; difference 14.4 percentage points","95% CI 7.5 to 21.3; P<0.001","40305708","10.1056/NEJMoa2413258",AR,"ESSENCE part 1","yes")
row(*m,"NCT03357380","1","67","48 (72 total)","NAFLD, MRI-PDFF >=10% and MRE 2.50-4.63 kPa","pdff","improved_unadjusted","estimated treatment ratio for liver fat 0.47 at week 48 (0.70 at week 24, 0.50 at week 72)","95% CI 0.36 to 0.60; P<0.0001","34570916","10.1111/apt.16608",AR,"semaglutide 0.4 mg daily; secondary endpoint; primary (MRE stiffness) not met","yes")
# ---- GLP-1R/GIPR
m=("GLP-1R/GIPR dual agonist","tirzepatide","GLP1R;GIPR","activate","CHEMBL4297839")
row(*m,"NCT04166773","2","190","52","biopsy MASH, fibrosis F2-F3","mash_resolution","improved","15 mg 62% vs 10%; difference 53 percentage points","95% CI 37 to 69; P<0.001 (registry logistic OR 19.63, 5.73 to 67.25)","38856224","10.1056/NEJMoa2401943",AR,"SYNERGY-NASH primary endpoint; 5 and 10 mg also P<0.001","yes")
row(*m,"NCT04166773","2","190","52","biopsy MASH, fibrosis F2-F3","fibrosis_improvement","improved_unadjusted","15 mg 51% vs 30%; difference 21 percentage points","paper 95% CI 1 to 42; registry logistic OR 2.47 (0.95 to 6.40) p=0.063; 5 mg registry p=0.025","38856224","10.1056/NEJMoa2401943",AR,"key secondary; paper (placebo-pattern imputation) and registry logistic analysis disagree on whether CI excludes null","yes")
row(*m,"NCT04166773","2","190","52","biopsy MASH, fibrosis F2-F3","pdff","improved_unadjusted","15 mg LS mean difference in absolute MRI-PDFF -10.02 percentage points","95% CI -13.17 to -6.87; p<0.001",NV,NV,REG,"secondary endpoint taken from registry results; not in abstract","yes")
# ---- GLP-1R/GCGR
m=("GLP-1R/GCGR dual agonist","survodutide","GLP1R;GCGR","activate","CHEMBL5314776")
row(*m,"NCT04771273","2","293 dosed (registry 295)","48","biopsy MASH, fibrosis F1-F3","fibrosis_improvement","improved_unadjusted","6.0 mg 44.2% vs 21.5% (registry actual-maintenance set); abstract 34% vs 22%","registry OR 3.37 (95% CI 1.52 to 7.45), p=0.0028, labelled nominal","38847460","10.1056/NEJMoa2401755",AR,"secondary endpoint; registry states p-values are nominal","yes")
row(*m,"NCT04771273","2","293 dosed (registry 295)","48","biopsy MASH, fibrosis F1-F3","pdff","improved_unadjusted",">=30% relative MRI-PDFF reduction: abstract 6.0 mg 57% vs 14%","registry OR 16.09 (95% CI 6.69 to 38.73), p<0.0001, labelled nominal","38847460","10.1056/NEJMoa2401755",AR,"secondary endpoint; primary was MASH improvement (not resolution), so no mash_resolution row","yes")
# ---- triple
m=("GLP-1R/GIPR/GCGR triple agonist","retatrutide","GLP1R;GIPR;GCGR","activate","CHEMBL5095485")
row(*m,"NCT04881760","2","98 (MASLD substudy of 338-participant obesity trial)","24 (endpoint; 48 treatment)","obesity/overweight with MASLD and liver fat >=10%","pdff","improved","mean relative change in liver fat at 24 weeks: 12 mg -82.4% vs placebo +0.3%","P<0.001 vs placebo, all doses","38858523","10.1038/s41591-024-03018-2",AR,"substudy primary objective; parent trial registered for obesity","yes")
# ---- FGF21
row("FGF21 analogue","pegozafermin",NV,"activate","CHEMBL4594466","NCT04929483","2b","222","24","biopsy NASH, fibrosis F2-F3","fibrosis_improvement","improved","30 mg weekly 26% vs pooled placebo 7%; difference 19 percentage points (44 mg every 2 weeks: 27%, difference 20)","30 mg: 95% CI 5 to 32, P=0.009; 44 mg Q2W: 5 to 35, P=0.008","37356033","10.1056/NEJMoa2304286",AR,"ENLIVEN; ChEMBL has no mechanism record; abstract says FGF21 analogue, receptors not named","yes")
row("FGF21 analogue","pegozafermin",NV,"activate","CHEMBL4594466","NCT04929483","2b","222","24","biopsy NASH, fibrosis F2-F3","mash_resolution","improved","30 mg weekly 23% vs 2%; difference 21 percentage points (15 mg 37%, 44 mg Q2W 26%)","30 mg: 95% CI 9 to 33; P not given in abstract","37356033","10.1056/NEJMoa2304286",AR,"co-primary; coded improved from CI excluding 0","yes")
row("FGF21 analogue","efruxifermin",NV,"activate","CHEMBL4594389","NCT05039450","2b","181","36 (96 secondary)","biopsy MASH compensated cirrhosis (F4)","fibrosis_improvement","not_improved","50 mg 19% vs 13%; difference 4 percentage points (week 96: 29% vs 11%, difference 16)","week 36: 95% CI -10 to 18, P=0.52; week 96: 2 to 30 (secondary)","40341827","10.1056/NEJMoa2502242",AR,"SYMMETRY; primary endpoint at week 36 failed","no")
row("FGF21 analogue","efruxifermin",NV,"activate","CHEMBL4594389","NCT04767529","2b","128","24","biopsy NASH, fibrosis F2-F3","fibrosis_improvement","improved","50 mg 41% vs 20% (28 mg 39%)","50 mg RR 2.2 (95% CI 1.0 to 5.0) p=0.036; 28 mg RR 2.3 (1.1 to 4.8) p=0.025; full analysis set 50 mg p=0.123","37802088","10.1016/S2468-1253(23)00272-8",AR,"HARMONY primary endpoint in liver biopsy analysis set; 96-week report PMID 40818852: 50 mg 49% vs 19%, p=0.0030","no")
row("FGF21 analogue","efruxifermin",NV,"activate","CHEMBL4594389","NCT05039450","2b","31","12","T2D with MASH F1-F3 on stable GLP-1RA","pdff","improved_unadjusted","hepatic fat fraction relative change -65% vs -10%","P<0.0001","38447814","10.1016/j.cgh.2024.02.022",AR,"SYMMETRY cohort D, add-on to GLP-1RA; secondary endpoint (primary was safety); MRI method for HFF not named in abstract","yes")
# ---- PPAR
m=("pan-PPAR agonist","lanifibranor","PPARA;PPARD;PPARG","activate","CHEMBL4091374")
row(*m,"NCT03008070","2b","247","24","noncirrhotic highly active NASH","mash_resolution","no_test_reported","1200 mg 49% vs 22% (800 mg 39%)","no p-value in abstract","34670042","10.1056/NEJMoa2036205",AR,"NATIVE; secondary endpoint; primary (SAF-A >=2-point decrease) met at 1200 mg P=0.007; phase 3 NATiV3 (NCT04849728) has no readout","yes")
row(*m,"NCT03008070","2b","247","24","noncirrhotic highly active NASH","fibrosis_improvement","no_test_reported","1200 mg 48% vs 29% (800 mg 34%)","no p-value in abstract","34670042","10.1056/NEJMoa2036205",AR,"NATIVE; secondary endpoint","yes")
m=("PPAR-gamma agonist","pioglitazone","PPARG","activate","CHEMBL595")
row(*m,"NCT00994682","4","101","78","biopsy NASH with prediabetes or T2D","mash_resolution","improved_unadjusted","51% resolution; treatment difference 32 percentage points","95% CI 13 to 51; P<0.001","27322798","10.7326/M15-1774",AR,"secondary endpoint; 45 mg/d; single centre; ChEMBL mechanism stored on salt CHEMBL1715","yes")
row(*m,"NCT00994682","4","101","78","biopsy NASH with prediabetes or T2D","fibrosis_improvement","improved_unadjusted","fibrosis score treatment difference -0.5","95% CI -0.9 to 0.0; P=0.039","27322798","10.7326/M15-1774",AR,"mean fibrosis score, not a >=1-stage responder endpoint","yes")
row(*m,"NCT00994682","4","101","78","biopsy NASH with prediabetes or T2D","pdff","improved_unadjusted","hepatic triglyceride 19% to 7%; treatment difference -7 percentage points","95% CI -10 to -4; P<0.001","27322798","10.7326/M15-1774",AR,"measured by proton MR spectroscopy, not MRI-PDFF","yes")
row(*m,"NCT00063622","3","247 (80 pioglitazone, 83 placebo)","96","nondiabetic NASH","alt","improved_unadjusted","ALT reduced vs placebo",
    "P<0.001","20427778","10.1056/NEJMoa0907929",AR,"PIVENS secondary endpoint; 30 mg/d; primary composite not met for pioglitazone (34% vs 19%, P=0.04, threshold 0.025)","yes")
row(*m,"NCT00063622","3","247 (80 pioglitazone, 83 placebo)","96","nondiabetic NASH","fibrosis_improvement","not_improved","fibrosis score not improved","P=0.12","20427778","10.1056/NEJMoa0907929",AR,"PIVENS; fibrosis score, not >=1-stage responder","no")
m=("vitamin E (antioxidant)","vitamin E","none (no ChEMBL mechanism)","not applicable","CHEMBL47 (alpha-tocopherol; no mechanism)")
row(*m,"NCT00063622","3","247 (84 vitamin E, 83 placebo)","96","nondiabetic NASH","alt","improved_unadjusted","ALT reduced vs placebo","P<0.001","20427778","10.1056/NEJMoa0907929",AR,"PIVENS secondary; 800 IU/d; primary composite histology improved 43% vs 19% P=0.001","yes")
row(*m,"NCT00063622","3","247 (84 vitamin E, 83 placebo)","96","nondiabetic NASH","fibrosis_improvement","not_improved","fibrosis score not improved","P=0.24","20427778","10.1056/NEJMoa0907929",AR,"PIVENS; fibrosis score, not >=1-stage responder; NASH resolution not in abstract","yes")
m=("PPAR-alpha/delta agonist","elafibranor","PPARA;PPARD","activate","CHEMBL3707395")
row(*m,"NCT02704403","3","2157 (F2-F3 analysis set 1070)","72","NASH with fibrosis F2-F3","mash_resolution","not_improved","120 mg 138/717 (19.2%) vs 52/353 (14.7%); net difference 0.043","95% CI -0.003 to 0.090; p=0.0659 (two-sided alpha 0.01)",NV,NV,REG,"RESOLVE-IT interim histology; trial terminated","yes")
row(*m,"NCT02704403","3","2157 (F2-F3 analysis set 1070)","72","NASH with fibrosis F2-F3","fibrosis_improvement","no_test_reported","120 mg 176/717 (24.5%) vs 79/353 (22.4%)","no analysis posted",NV,NV,REG,"RESOLVE-IT secondary endpoint","yes")
row(*m,"NCT01694849","2b","276 (93/91/92; registry 275)","52","noncirrhotic NASH","alt","improved_unadjusted","liver enzymes significantly reduced with 120 mg vs placebo (no ALT value in abstract)",NV,"26874076","10.1053/j.gastro.2016.01.038",AR,"GOLDEN-505; primary NASH resolution endpoint not met in ITT","yes")
m=("PPAR-alpha/gamma agonist","saroglitazar","PPARA;PPARG","activate","CHEMBL4297530")
row(*m,"NCT03061721","2","106","16","NAFLD/NASH with ALT >=50 U/L and BMI >=25","alt","improved","4 mg LS mean % change in ALT -45.8% vs +3.4% placebo","P<0.001","33811367","10.1002/hep.31843",AR,"EVIDENCES IV primary endpoint; 1 and 2 mg also P<0.001","yes")
row(*m,"NCT03061721","2","106","16","NAFLD/NASH with ALT >=50 U/L and BMI >=25","pdff","improved_unadjusted","4 mg liver fat content change -19.7% vs +4.1% placebo","P<0.05","33811367","10.1002/hep.31843",AR,"secondary endpoint, MRI-PDFF","yes")
# ---- FXR
m=("FXR agonist","obeticholic acid","NR1H4","activate","CHEMBL566315")
row(*m,"NCT02548351","3","931 F2-F3 primary analysis (registry 2477)","78","NASH, fibrosis F2-F3","fibrosis_improvement","improved","25 mg 23% vs 12% (10 mg 18%)","25 mg p=0.0002; 10 mg p=0.045","31813633","10.1016/S0140-6736(19)33041-7",AR,"REGENERATE month-18 interim; 2023 consensus re-read (PMID 37517454) 22.4% vs 9.6%, p<0.0001","yes")
row(*m,"NCT02548351","3","931 F2-F3 primary analysis (registry 2477)","78","NASH, fibrosis F2-F3","mash_resolution","not_improved","25 mg 12% vs 8%","p=0.13","31813633","10.1016/S0140-6736(19)33041-7",AR,"REGENERATE; 2023 re-read 6.5% vs 3.5%, p=0.093","yes")
m=("FXR agonist","cilofexor","NR1H4","activate","CHEMBL4297613")
row(*m,"NCT03449446","2b","392","48","NASH with bridging fibrosis or compensated cirrhosis (F3-F4)","fibrosis_improvement","not_improved","cilofexor 30 mg 12% vs placebo 11%","P=0.96","33169409","10.1002/hep.31622","PubMed abstract; registry record (NCT matched by title and design, not stated in abstract)","ATLAS monotherapy arm","no")
row(*m,"NCT02854605","2","140","24","noncirrhotic NASH by MRI-PDFF >=8% and stiffness","pdff","improved_unadjusted","100 mg median relative MRI-PDFF change -22.7% vs +1.9%","P=0.003 (30 mg P=0.17)","32115759","10.1002/hep.31205",AR,"abstract does not label this the primary endpoint","no")
m=("FXR agonist","tropifexor","NR1H4","activate","CHEMBL4298169")
row(*m,"NCT02855164","2","350 (Parts A+B 198, Part C 152)","12 (endpoint; 48 Part C)","NASH","alt","no_test_reported","week 12 ALT change: 200 ug -23.0 U/L vs placebo -8.3 U/L","no p-value in abstract","36797481","10.1038/s41591-022-02200-8",AR,"FLIGHT-FXR co-primary dose-response; abstract states decreases were greater than placebo","yes")
row(*m,"NCT02855164","2","350 (Parts A+B 198, Part C 152)","12 (endpoint; 48 Part C)","NASH","pdff","no_test_reported","week 12 hepatic fat fraction change: 200 ug -39.41% vs placebo -10.77%","no p-value in abstract","36797481","10.1038/s41591-022-02200-8",AR,"HFF by MRI; method name not given in abstract","yes")
# ---- ACC
m=("ACC inhibitor","firsocostat","ACACA;ACACB","inhibit","CHEMBL3407547")
row(*m,"NCT03449446","2b","392","48","NASH with bridging fibrosis or compensated cirrhosis (F3-F4)","fibrosis_improvement","not_improved","firsocostat 20 mg 12% vs placebo 11%","P=0.94","33169409","10.1002/hep.31622","PubMed abstract; registry record (NCT matched by title and design, not stated in abstract)","ATLAS monotherapy arm; ChEMBL action type negative allosteric modulator","yes")
row(*m,"NCT02856555","2","126 (registry 127)","12","NASH by MRI-PDFF >=8% and stiffness","pdff","improved_unadjusted","20 mg: >=30% PDFF response 48% vs 15%; median relative change -29% vs -8%","P=0.004 (responders); P=0.002 (median)","30059671","10.1053/j.gastro.2018.07.027",AR,"GS-0976; abstract does not name a primary endpoint","yes")
# ---- DGAT2
row("DGAT2 inhibitor","ervogastat","DGAT2","inhibit","CHEMBL4760665","NCT04321031","2","255","48","biopsy MASH, fibrosis F2-F3","fibrosis_improvement","not_improved","no ervogastat monotherapy group exceeded placebo (no numbers in abstract)","composite primary not met for monotherapy (e.g. 75 mg 52% vs 38%, difference 0.14, 90% CI -0.04 to 0.32)","40753985","10.1016/S2468-1253(25)00128-1",AR,"MIRNA; ChEMBL has no mechanism; DGAT2 target from PMID 34635855 (PF-06865571 = ervogastat per ChEMBL synonym)","yes")
row("DGAT2 inhibitor","ervogastat","DGAT2","inhibit","CHEMBL4760665","NCT04321031","2","255","48","biopsy MASH, fibrosis F2-F3","mash_resolution","no_test_reported","abstract: all experimental groups showed greater effects than placebo; no numbers",NV,"40753985","10.1016/S2468-1253(25)00128-1",AR,"MIRNA secondary component","yes")
row("DGAT2 inhibitor","IONIS-DGAT2Rx","DGAT2","inhibit (antisense knockdown)",NV,"NCT03334214","2","44","13","T2D with liver fat >=10%","pdff","improved","absolute liver fat change -5.2% vs -0.6%; difference -4.2 percentage points","95% CI -7.8 to -0.5; p=0.026","32553151","10.1016/S2468-1253(20)30186-2",AR,"primary pharmacodynamic endpoint, per-protocol; ION224 phase 2 (PMID 40849139, n=160) met a NAS-based primary not among the four endpoints","yes")
# ---- SCD1
m=("SCD1 partial inhibitor","aramchol","SCD","inhibit","CHEMBL6068507 (aramchol meglumine; no mechanism)")
p=("NCT02279524","2b","247","52","NASH")
row(*m,*p,"pdff","not_improved","600 mg placebo-corrected liver triglyceride change -3.1","95% CI -6.4 to 0.2; P=0.066","34621052","10.1038/s41591-021-01495-3",AR,"ARREST primary; measured by MR spectroscopy, not MRI-PDFF; target SCD from abstract text","yes")
row(*m,*p,"mash_resolution","not_improved","600 mg 16.7% (13/78) vs 5% (2/40); OR 4.74","95% CI 0.99 to 22.7; formal testing precluded","34621052","10.1038/s41591-021-01495-3",AR,"key secondary after failed primary","yes")
row(*m,*p,"fibrosis_improvement","not_improved","600 mg 29.5% vs 17.5%; OR 1.88","95% CI 0.7 to 5.0","34621052","10.1038/s41591-021-01495-3",AR,"key secondary after failed primary","yes")
row(*m,*p,"alt","improved_unadjusted","600 mg placebo-corrected ALT change -29.1 IU/L","95% CI -41.6 to -16.5","34621052","10.1038/s41591-021-01495-3",AR,"key secondary; hierarchy had stopped","yes")
# ---- FASN
m=("FASN inhibitor","denifanstat","FASN","inhibit","CHEMBL3661754")
row(*m,"NCT04906421","2b","168","52","biopsy MASH, fibrosis F2-F3","mash_resolution","improved","26% vs 11% (MASH resolution with >=2-point NAS improvement); risk difference 13.0%","95% CI 0.7 to 25.3; p=0.0173","39396529","10.1016/S2468-1253(24)00246-2",AR,"FASCINATE-2 primary; ChEMBL has no mechanism; FASN from PMID 34310978 title","yes")
row(*m,"NCT04906421","2b","168","52","biopsy MASH, fibrosis F2-F3","fibrosis_improvement","improved_unadjusted",NV,NV,"39396529","10.1016/S2468-1253(24)00246-2",AR,"abstract interpretation states statistically significant fibrosis improvement but gives no numbers","yes")
# ---- ASK1
m=("ASK1 inhibitor","selonsertib","MAP3K5","inhibit","CHEMBL3916717")
row(*m,"NCT03053063","3","877 (354/351/172; registry 883)","48","NASH compensated cirrhosis (F4)","fibrosis_improvement","not_improved","18 mg 14% vs 13%","p=0.56","32147362","10.1016/j.jhep.2020.02.027",AR,"STELLAR-4 primary","yes")
row(*m,"NCT03053050","3","802 (322/321/159; registry 808)","48","NASH bridging fibrosis (F3)","fibrosis_improvement","not_improved","18 mg 10% vs 13%","p=0.49","32147362","10.1016/j.jhep.2020.02.027",AR,"STELLAR-3 primary","no")
row(*m,"NCT03053063","3","877 (354/351/172; registry 883)","48","NASH compensated cirrhosis (F4)","alt","not_improved","no significant effect on liver biochemistry (no value in abstract)",NV,"32147362","10.1016/j.jhep.2020.02.027",AR,"statement covers both STELLAR trials","yes")
# ---- caspase
m=("pan-caspase inhibitor","emricasan","CASP1;CASP2;CASP3;CASP4;CASP5;CASP6;CASP7;CASP8;CASP9;CASP10;CASP12;CASP14","inhibit","CHEMBL197672")
row(*m,"NCT02686762","2","318","72","NASH, fibrosis F1-F3","fibrosis_improvement","not_improved","50 mg 12.3% vs placebo 19.0% (5 mg 11.2%); OR 0.588","p=0.972","31887369","10.1016/j.jhep.2019.11.024",AR,"ENCORE-NF primary; authors state emricasan may have worsened fibrosis and ballooning (not a significant test)","yes")
row(*m,"NCT02686762","2","318","72","NASH, fibrosis F1-F3","mash_resolution","not_improved","50 mg 6.6% vs 10.5% (5 mg 3.7%); OR 0.613","p=0.335","31887369","10.1016/j.jhep.2019.11.024",AR,"ENCORE-NF secondary","yes")
row(*m,"NCT02077374","2","38","4","NAFLD with ALT >=1.5x ULN","alt","improved","ALT decreased vs placebo (no value in abstract)","day 28 P=0.02; day 7 P<0.0001","30430605","10.1111/apt.15030",AR,"primary endpoint; 25 mg twice daily; ENCORE-NF ALT not in abstract","yes")
# ---- CCR2/5
m=("CCR2/CCR5 antagonist","cenicriviroc","CCR2;CCR5","inhibit","CHEMBL2110727")
row(*m,"NCT03028740","3","1778 (part 1: 1293)","52","NASH, fibrosis F2-F3","fibrosis_improvement","not_improved","22.3% vs 25.5%; OR 0.84","95% CI 0.63 to 1.10; P=0.21","37061109","10.1016/j.cgh.2023.04.003",AR,"AURORA part 1 primary; trial terminated","yes")
row(*m,"NCT03028740","3","1778 (part 1: 1293)","52","NASH, fibrosis F2-F3","mash_resolution","not_improved","23.0% vs 27.2%","P=0.21","37061109","10.1016/j.cgh.2023.04.003",AR,"AURORA part 1","yes")
# ---- galectin-3
row("galectin-3 inhibitor","belapectin","LGALS3","inhibit","CHEMBL4297577","NCT02462967","2b","162","52","NASH cirrhosis with portal hypertension (HVPG >=6)","fibrosis_improvement","not_improved","no significant effect on fibrosis (no value in abstract)",NV,"31812510","10.1053/j.gastro.2019.11.296",AR,"secondary histology; primary HVPG not met; NAVIGATE (PMID 42065864, n=357) reported varices only","yes")
# ---- LOXL2
row("LOXL2 antibody","simtuzumab","LOXL2","inhibit","CHEMBL2109667","NCT01672866","2b","219 (registry 222)","96","NASH bridging fibrosis (F3)","fibrosis_improvement","not_improved","fibrosis stage not significantly decreased; primary collagen change 125 mg vs placebo -0.4%","primary: 95% CI -1.5 to 0.8, P=0.52","29990488","10.1053/j.gastro.2018.07.006",AR,"primary endpoint was morphometric collagen, not >=1-stage; cirrhosis trial NCT01672879 (n=258) primary HVPG also null","yes")
# ---- SGLT2
m=("SGLT2 inhibitor","dapagliflozin","SLC5A2","inhibit","CHEMBL429910")
row(*m,"NCT03723252","3","154","48","biopsy MASH with or without T2D (China)","mash_resolution","improved_unadjusted","23% (18/78) vs 8% (6/76); RR 2.91","95% CI 1.22 to 6.97; P=0.01","40467095","10.1136/bmj-2024-083735",AR,"DEAN secondary; primary (MASH improvement) met P=0.006","yes")
row(*m,"NCT03723252","3","154","48","biopsy MASH with or without T2D (China)","fibrosis_improvement","improved_unadjusted","45% (35/78) vs 20% (15/76); RR 2.25","95% CI 1.35 to 3.75; P=0.001","40467095","10.1136/bmj-2024-083735",AR,"DEAN secondary","yes")
m=("SGLT2 inhibitor","empagliflozin","SLC5A2","inhibit","CHEMBL2107830")
row(*m,"NCT04642261","4","98","52","nondiabetic MASLD, MRI-PDFF >=5%","pdff","improved","median MRI-PDFF change -2.49% vs -1.43%","p=0.025","38536017","10.1097/HEP.0000000000000855",AR,"primary endpoint","yes")
row(*m,"NCT04642261","4","98","52","nondiabetic MASLD, MRI-PDFF >=5%","alt","not_improved","ALT drop >=17 U/L: 16.3% vs 12.2%","p=0.564","38536017","10.1097/HEP.0000000000000855",AR,"secondary responder endpoint, not mean ALT change","yes")
# ---- MSDC-0602K
m=("insulin sensitizer (MPC modulator per candidate list)","MSDC-0602K",NV,"inhibit (not verified)","CHEMBL4593380 (no mechanism)")
p=("NCT02784444","2b","392","52","biopsy NASH, fibrosis F1-F3")
row(*m,*p,"mash_resolution","not_improved","no significant effect (no value in abstract)",NV,"31697972","10.1016/j.jhep.2019.10.023",AR,"EMMINENCE secondary; abstract names no molecular target; primary (NAS-based) not met, 250 mg OR 1.64 (0.83 to 3.27)","yes")
row(*m,*p,"fibrosis_improvement","not_improved","no significant effect (no value in abstract)",NV,"31697972","10.1016/j.jhep.2019.10.023",AR,"EMMINENCE secondary","yes")
row(*m,*p,"alt","improved_unadjusted","two highest doses significantly reduced liver enzymes (no value in abstract)",NV,"31697972","10.1016/j.jhep.2019.10.023",AR,"exploratory endpoint","yes")
# ---- KHK
row("KHK inhibitor","PF-06835919","KHK","inhibit","CHEMBL4549658","NCT03969719","2","164","16","NAFLD (MRI-PDFF >=8%) with T2D on metformin","pdff","improved","300 mg LS mean % change in whole liver fat -19.13% vs placebo -5.26%","90% CI -25.51 to -12.20 (drug arm); P=0.0288 vs placebo","36515213","10.1111/dom.14946","PubMed abstract; registry record (NCT matched by drug, n=164, 16 weeks; not stated in abstract)","primary endpoint; ChEMBL has no mechanism; KHK from abstract text and ChEMBL synonym 'ketohexokinase inhibitor 1'","yes")
# ---- FGF19
row("FGF19 analogue","aldafermin","FGF1;FGFR4;KLB","activate","CHEMBL4298055","NCT03912532","2b","171","24","biopsy NASH, fibrosis F2-F3","fibrosis_improvement","not_improved","3 mg 30% vs 19%; difference 10 percentage points","90% CI -9 to 30, p=0.12; dose-response p=0.55","35325622","10.1016/S2468-1253(22)00017-6",AR,"ALPINE 2/3 primary; ChEMBL target CHEMBL4630720 lists FGF1 (P05230) with KLB, copied as listed","yes")
# ---- A3AR
row("adenosine A3 receptor agonist","namodenoson","ADORA3","activate","CHEMBL431733","NCT02927314","2","60","12","NAFLD with ALT >=60 IU/L","alt","not_improved","25 mg twice daily ALT change vs placebo (no value in abstract)","P=0.066","34671996","10.1111/apt.16664",AR,"main efficacy endpoint; AST P=0.03","yes")
# ---- AOC3
row("AOC3 (VAP-1) inhibitor","BI 1467335","AOC3","inhibit","CHEMBL3910400 (no mechanism)","NCT03166735","2","114","12","clinical evidence of NASH","alt","improved_unadjusted","ALT as % of baseline: 10 mg 77.57 vs placebo 92.66 (means)","MCPMod dose-response p=0.0212 (one-sided alpha 0.05)",NV,NV,REG,"secondary endpoint; target from registry description ('AOC3 inhibition')","yes")
# ---- PNPLA3 knockdown
m=("PNPLA3 antisense knockdown","AZD2693","PNPLA3","inhibit (antisense knockdown)",NV)
row(*m,"NCT05809934","2b","216 started (73/69/74; registry 220)","52","non-cirrhotic MASH, PNPLA3 148M carriers","mash_resolution","not_improved","50 mg 1/69 vs placebo 5/74; OR 0.20 (25 mg 4/73, OR 0.83)","50 mg 95% CI 0.02 to 1.78, p=0.116; 25 mg 0.21 to 3.23, p=0.783",NV,NV,REG,"FORTUNA primary; no publication found; ALT >3x ULN in 33/69 (50 mg) vs 10/74 placebo (safety table, not tested)","yes")
row(*m,"NCT05809934","2b","216 started (73/69/74; registry 220)","52","non-cirrhotic MASH, PNPLA3 148M carriers","fibrosis_improvement","not_improved","50 mg 6/69 vs placebo 8/74; OR 0.80 (25 mg 12/73, OR 1.64)","50 mg 95% CI 0.26 to 2.44, p=0.693; 25 mg 0.62 to 4.37, p=0.320",NV,NV,REG,"FORTUNA secondary","yes")
row(*m,"NCT04483947","1","not verified (registry enrollment 74)","12","MRI-PDFF >=7%, PNPLA3 148M homozygous","pdff","no_test_reported","placebo-corrected change in hepatic steatosis at week 12: 50 mg -12.2%, 25 mg -7.6%","no p-value in abstract","39798707","10.1016/j.jhep.2024.12.046",AR,"phase 1 multiple ascending dose; units (relative vs absolute) not stated in abstract","yes")

with open(os.path.join(D,"trial_outcomes.tsv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(R)

EX=[("HSD17B13 knockdown","rapirosiran (ALN-HSD); GSK4532990; ARO-HSD; AZD7503; INI-822","HSD17B13","phase 1 to 2a only; no placebo-controlled readout on the four endpoints",
     "Rapirosiran phase 1 (NCT04565717, PMID 40581300) reported safety and mRNA knockdown only; GSK4532990 SKYLINE (NCT06104319) phase 2, no placebo arm (sequential, open label), mRNA/protein endpoints; GSK4532990 phase 2b HORIZON (NCT05583344) active, not recruiting, no results; ARO-HSD NCT04202354 phase 1/2 results are PK only in the outcomes read"),
    ("CIDEB knockdown","ALN-CIDEB","CIDEB","no randomized readout verified",
     "NCT06836609 phase 1/2 (n=118) active, no posted results; NEJM correspondence PMID 43001761 (2026) has no abstract and was not read; phase 2 NCT07631637 and phase 3 NCT07854782 recruiting"),
    ("MTARC1 knockdown","none found","MTARC1","no interventional drug trial found","ClinicalTrials.gov searches for MTARC1, ARO-MARC1 and mitochondrial amidoxime returned no drug trial"),
    ("INHBE knockdown","ARO-INHBE","INHBE","no MASLD readout","NCT06700538 phase 1/2 in obesity (active, not recruiting); NCT07859813 phase 2 in MASLD not yet recruiting"),
    ("GPAM knockdown","none found","GPAM","no interventional drug trial found","ClinicalTrials.gov searches for GPAM and ALN-GPAM returned no drug trial"),
    ("PNPLA3 siRNA (other programs)","JNJ-75220795; LY3849891; ALN-PNP","PNPLA3","phase 1 or recruiting","JNJ-75220795 phase 1 (NCT04844450); LY3849891 phase 1 (NCT05395481); ALN-PNP phase 2a recruiting (NCT07527910). PNPLA3 mechanism is represented by AZD2693 FORTUNA"),
    ("glucocorticoid receptor modulator","miricorilant","NR3C1","no randomized controlled readout","Phase 1b (PMID 42929711) had no placebo arm in abstract; phase 2b MONARCH NCT06108219 not read for results"),
    ("MOGAT2 inhibitor","BMS-963272","MOGAT2","not searched to completion","Listed in candidate file as phase 1-2; no MASLD randomized readout found in the time available"),
    ("RORA","candidate list only","RORA","phase 1 per candidate list; not searched","no MASLD trial opened"),
    ("CB1 (CNR1)","candidate list only","CNR1","not searched","no MASLD trial opened"),
    ("GLP-1R/GCGR (other drugs)","pemvidutide; efinopegdutide; cotadutide","GLP1R;GCGR","mechanism already represented","survodutide phase 2 chosen by the fixed rule; other drugs not extracted"),
    ("THR-beta (other drug)","VK2809","THRB","mechanism already represented","resmetirom phase 3 chosen by the fixed rule"),
    ("ACC inhibitor (other drug)","clesacostat (PF-05221304)","ACACA;ACACB","mechanism already represented","phase 2a PDFF trial PMID 34635855 is lower phase than firsocostat phase 2 (NCT02856555)"),
    ("DGAT2 antisense (other)","ION224","DGAT2","no row on the four endpoints","ION224-CS2 (NCT04932512, PMID 40849139) primary was NAS-based histologic improvement, met at 90 and 120 mg; resolution/fibrosis/PDFF values not in abstract"),
    ("saroglitazar NASH histology trial","saroglitazar","PPARA;PPARG","no abstract","PMID 33152542 is a research letter without abstract; not read"),
]
with open(os.path.join(D,"excluded.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["mechanism","drug","target_genes","status","reason_and_sources"]); w.writerows(EX)
print(len(R),"rows;",len(EX),"excluded")
