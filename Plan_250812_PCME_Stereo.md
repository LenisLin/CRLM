# **Peritumor Cholangiocyte Microenvironment (PCME) Analysis in CRLM Recurrence**
## **Executive Project Plan for Spatial Transcriptomics Validation**

---

## **1. Background & Hypothesis**

### **1.1 Research Context**
Colorectal liver metastases (CRLM) recurrence prediction remains challenging. Our previous IMC analysis of 35 patients revealed that **Peritumor Cholangiocyte Microenvironment (PCME)** serves as a novel "soil" biomarker for recurrence prediction.

### **1.2 Key IMC Discoveries**
- **PCME-I (Immune-infiltrated)**: High immune infiltration → **Non-early relapse**
- **PCME-S (Stromal-fibrotic)**: High Vimentin+/Collagen-I+ → **Early relapse**
- **Clinical Performance**: AUC=0.75, log-rank p=0.0016, C-index=0.8

### **1.3 Central Hypothesis**
**Cholangiocyte stress reprogramming and their spatial microenvironment composition determine CRLM recurrence risk through distinct immune vs stromal "soil" states.**

### **1.4 Spatial Transcriptomics Objective**
Validate and extend PCME findings using spatial transcriptomics to discover:
- Molecular programs driving cholangiocyte reprogramming
- Spatial organization principles of PCME-I vs PCME-S
- Therapeutic targets for microenvironment modulation

---

## **2. Dataset Description**

### **2.1 Spatial Transcriptomics Cohort**
- **Total samples**: 9 FDZS slides (PCME-enriched regions)
- **Early recurrence (Ax)**: 5 samples (A1, A2, A4, A5, A6, A7)
- **Non-early recurrence (Cx)**: 3 samples (C5, C6, C7)
- **Technology**: Stereo-seq, bin50 resolution (multiple cells resolution)

### **2.2 Cell Type Annotation**
- **Deconvolution method**: Cell2location using integrated scRNA-seq reference
- **Cell types**: 44 subtypes across major lineages:
  - **Epithelial**: General epithelial, Hepatocytes
  - **Immune**: T cells (CD4T, CD8T, Tregs), B cells, Macrophages, NK cells, DCs
  - **Stromal**: CAFs, Fibroblasts, Endothelial cells
  - **Tumor**: 5 malignant cell programs (EMT, Glycolysis, etc.)

---

## **3. Detailed Analysis Framework**

### **Phase 1: Sample-Level Characterization & Hypothesis Validation**

#### **Step 1.1: Data Description and Quality Control**
```
Objective: Comprehensive data overview and quality assessment
Methods:
- Load 9 FDZS slides with cell2location abundance data
- Sample-level statistics: total cells, cell type diversity
- Visualization: Sample composition overview, spatial maps

Expected Output:
- Cell type composition summary (A vs C groups)
- Spatial distribution maps for each sample
```

#### **Step 1.2: Global Cell Type Differential Analysis**
```
Objective: Compare cellular landscape and test H1 - Ax samples have higher stromal abundance
Methods:
- Statistical comparison: A-group vs C-group for all 44 cell types
- Multiple testing correction (FDR < 0.05, BH adjustment, wilcox test)
- Effect size calculation (fold-change)
- Visualization: Stacked bar plots, violin plots, heatmaps

Expected Output:
- Statistical summary of cell type differences (A vs C)
- Ranked list of recurrence-associated cell types
- Validation of IMC-identified populations (Vimentin+, Collagen-I+)
```

#### **Step 1.3: Stromal vs Immune Balance Analysis**
```
Objective: Test stromal/immune hypothesis with specific metrics
Methods:
- Define stromal signature: CAF_*, Fibro_*, Collagen markers
- Define immune signature: CD8T_*, CD4T_*, Macro_* (M1-like)*
- Calculate Stromal/Immune ratios per sample
- Statistical testing: Wilcoxon rank-sum test (A vs C)

Key Metrics:
- Stromal Index = (CAF + Fibroblast abundance) / Total cells
- Immune Index = (T cells + M1 macrophages) / Total cells
- S/I Ratio = Stromal Index / Immune Index

Expected Output:
- Validation of stromal enrichment in A-group samples
- Quantitative S/I ratio differences between groups
- Statistical significance testing results
```

---

### **Phase 2: Cholangiocyte Identification & Characterization**

#### **Step 2.1: Cholangiocyte Spot Identification**
```
Objective: Define cholangiocyte-enriched spots in spatial data
Methods:
- Cholangiocyte signature genes: CK7, EPCAM, GGT, AQP
- Map signature scores to "Epithelial" cell abundance spots
- Threshold definition: Top 20% epithelial spots with cholangiocyte markers
- Validate with known cholangiocyte anatomical locations

Expected Output:
- Cholangiocyte-enriched spot identification across all samples
- Signature score validation and threshold determination
- Spatial maps highlighting cholangiocyte regions
```

#### **Step 2.2: CAIX Expression Analysis in Cholangiocytes**
```
Objective: Test H3 - CAIX+ cholangiocytes associate with recurrence
Methods:
- Map CA9 (CAIX) gene expression to cholangiocyte spots
- Compare CAIX levels: A-group vs C-group cholangiocytes
- Identify CAIX+ vs CAIX- cholangiocyte subpopulations

Key Analysis:
- CAIX expression distribution in cholangiocyte spots
- Statistical testing: A vs C group CAIX levels
- CAIX+ spot percentage per sample group

Expected Output:
- Validation of CAIX upregulation in A-group cholangiocytes
- CAIX+ cholangiocyte spatial distribution patterns
- Correlation between CAIX and recurrence status
```

#### **Step 2.3: Cholangiocyte Differential Gene Expression**
```
Objective: Identify molecular programs in recurrence-associated cholangiocytes
Methods:
- Extract gene expression from cholangiocyte-enriched spots
- DEG analysis: A-group vs C-group cholangiocytes
- Pathway enrichment analysis of DEGs
- Focus on stress response, hypoxia, inflammatory pathways

Expected Output:
- List of differentially expressed genes in A vs C cholangiocytes
- Pathway enrichment results (stress, hypoxia, inflammation)
- Molecular signatures of recurrence-associated cholangiocyte reprogramming
```

---

### **Phase 3: PCME Spatial Organization Analysis**

#### **Step 3.1: Spatial Distance Analysis Framework**
```
Objective: Analyze cell population distributions around cholangiocytes
Methods (adapted from attached PDF methodology):
- Define cholangiocyte spots as spatial anchors (0μm)
- Calculate distance bins: 0-50μm, 50-100μm, 100-150μm, 150-200μm
- Map cell type abundances at each distance interval
- Compare distance profiles: A-group vs C-group

Distance Analysis Protocol:
1. Identify cholangiocyte anchor spots per sample
2. Calculate spot-to-spot distances for all other spots
3. Assign spots to distance bins relative to nearest cholangiocyte
4. Average cell type abundances within each distance bin
5. Generate distance-abundance curves for each cell type

Expected Output:
- Distance-dependent cell abundance profiles
- Spatial organization differences between A and C groups
- Identification of critical distance ranges for recurrence association
```

#### **Step 3.2: PCME-I vs PCME-S Spatial Signatures**
```
Objective: Define and validate spatial microenvironment subtypes
Methods:
- Identify spots within 100μm of cholangiocytes (PCME regions)
- Calculate immune and stromal enrichment scores per PCME region
- Classify PCME regions: PCME-I (immune-rich) vs PCME-S (stromal-rich)
- Compare spatial organization patterns

PCME Classification:
- PCME-I: High CD8T + CD4T + M1-Macrophage abundance, Low CAF + Fibro abundance
- PCME-S: High CAF + Fibro + Collagen signature, Low immune abundance
- Boundary definition: Median split or clustering-based classification

Expected Output:
- PCME-I vs PCME-S classification for all cholangiocyte neighborhoods
- Validation that A-group enriches for PCME-S, C-group enriches for PCME-I
- Spatial maps showing PCME subtype distributions
```

---

### **Phase 4: Molecular Mechanism Discovery**

#### **Step 4.1: PCME-Specific Pathway Analysis**
```
Objective: Identify molecular pathways driving PCME-I vs PCME-S formation
Methods:
- Extract gene expression from PCME-I and PCME-S regions
- Differential expression analysis: PCME-S vs PCME-I
- Pathway enrichment analysis (GO, KEGG, Reactome, Hallmarks)
- Focus on ECM remodeling, immune regulation, angiogenesis pathways

Expected Pathway Categories:
- PCME-S enriched: ECM organization, TGF-β signaling, collagen synthesis
- PCME-I enriched: T cell activation, interferon signaling, antigen presentation

Expected Output:
- Comprehensive pathway enrichment results for PCME subtypes
- Molecular signatures defining PCME-I vs PCME-S states
- Potential therapeutic targets for microenvironment modulation
```

#### **Step 4.2: Recurrence Prediction Signature Development**
```
Objective: Develop spatial transcriptomic signature for recurrence prediction
Methods:
- Integrate cholangiocyte stress markers + PCME spatial composition
- Feature selection: Top predictive genes and spatial metrics
- Cross-validation within 9-sample cohort
- Compare with IMC-derived PCME signatures

Signature Components:
- CAIX expression in cholangiocytes
- Stromal/Immune ratio within 100μm of cholangiocytes  
- Specific pathway activity scores (ECM vs immune)
- Spatial organization metrics

Expected Output:
- Integrated spatial transcriptomic recurrence signature
- Performance metrics and validation results
- Comparison with existing IMC-based PCME classification
```

#### **Step 4.3: Therapeutic Target Identification**
```
Objective: Identify druggable targets for PCME modulation
Methods:
- Druggable gene analysis in PCME-S vs PCME-I DEGs
- Ligand-receptor interaction analysis
- Connection to known drug databases (DrugBank, GDSC)
- Prioritize targets based on clinical tractability

Target Categories:
- ECM remodeling inhibitors (for PCME-S intervention)
- Immune activation enhancers (for PCME conversion)
- Cholangiocyte stress response modulators

Expected Output:
- Ranked list of therapeutic targets with druggability scores
- Proposed combination therapy strategies
- Mechanistic hypotheses for PCME-targeted interventions
```

---

## **4. Expected Deliverables**

### **4.1 Scientific Outputs**
1. **Validation Report**: Confirmation of IMC findings at transcriptional level
2. **Spatial Atlas**: Comprehensive maps of PCME organization patterns  
3. **Molecular Signatures**: Gene expression programs defining PCME subtypes
4. **Pathway Networks**: Regulatory circuits controlling PCME formation
5. **Therapeutic Targets**: Druggable pathway components for clinical translation

### **4.2 Technical Innovations**
1. **Spatial Analysis Pipeline**: Cholangiocyte-centered distance analysis framework
2. **PCME Classification**: Spatial transcriptomic method for microenvironment typing
3. **Integration Strategy**: Multi-modal validation approach (IMC + ST)

### **4.3 Clinical Translation Potential**
1. **Biomarker Validation**: Spatial transcriptomic confirmation of PCME prognostic value
2. **Mechanistic Insights**: Understanding of how normal liver tissue contributes to metastatic soil
3. **Therapeutic Strategy**: Rational basis for microenvironment-targeted therapy

---

## **5. Success Criteria**

### **5.1 Primary Validation Goals**
- [ ] Confirm higher stromal abundance in A-group samples (p < 0.05)
- [ ] Validate CAIX upregulation in A-group cholangiocytes (p < 0.05)  
- [ ] Demonstrate spatial organization differences between A and C groups
- [ ] Identify molecular pathways distinguishing PCME-I vs PCME-S

### **5.2 Discovery Goals**
- [ ] Define cholangiocyte reprogramming molecular programs
- [ ] Map therapeutic targets for PCME modulation
- [ ] Develop spatial transcriptomic recurrence prediction signature
- [ ] Establish framework for normal tissue "soil" analysis

---