# IMC Spatial Analysis Pipeline for CRLM Recurrence Prediction

## 1. Data Preparation and Quality Control

### 1.1 Spatial Graph Construction
- [ ] **Multiple spatial graph generation**
  - kNN graphs (k = 10, 20, 40 for different scales)
  - Expansion graphs (threshold = 4μm, 10μm, 20μm)
  - Delaunay triangulation graphs
  - steinbock neighborhood graphs
- [ ] **Parameter optimization**
  - Graph construction parameter sweep
  - Tissue density-adjusted parameters
  - Region-specific optimization (CT vs IM vs TAT)

### 1.2 Spatial Visualization Setup
- [ ] **Basic spatial plotting**
  - Cell centroid visualization with tissue context
  - Spatial graph edge visualization
  - Multi-channel marker expression overlay
- [ ] **Quality assessment**
  - Cell density distribution across images
  - Spatial graph connectivity metrics
  - Regional boundary definition validation

---

## 2. Individual Cell Type Interaction Analysis

### 2.1 Permutation-Based Interaction Testing

#### 2.1.1 Standard Interaction Analysis
- [ ] **Overall interaction patterns**
  - `testInteractions()` across all samples
  - Interaction vs avoidance quantification
  - Statistical significance testing (p-value correction)
- [ ] **Clinical stratification**
  - Relapse vs non-relapse interaction comparison
  - Treatment-specific interaction patterns
  - Regional interaction differences (CT vs IM vs TAT)

### 2.2 Malignant Cell-Microenvironment Interactions

#### 2.2.1 Treatment-Specific Analysis
- [ ] **EC_GLUT1 interaction profiling**
  - Chemotherapy-only cohort interactions
  - Combination therapy cohort interactions
  - Treatment-dependent interaction rewiring
- [ ] **Other malignant subtypes**
  - TC_CAIX, TC_EpCAM, TC_Ki67, TC_VEGF interaction patterns
  - Regional interaction specificity (TC vs IM)
  - Prognostic interaction signatures

#### 2.2.2 Invasion Biology Focus
- [ ] **Tumor margin interaction analysis**
  - EMT-activated cell interactions in IM
  - Tumor budding cell neighborhood characterization
  - Invasion-competence spatial signatures

---

## 3. Cellular Neighborhood (CN) Analysis

### 3.1 CN Detection and Optimization

#### 3.1.1 Multi-Method CN Detection
- [ ] **Cell type-based aggregation**
  - `aggregateNeighbors()` with celltype fractions
  - k-means clustering (k = 6, 10, 15, 20)
  - Parameter sweep and stability analysis

#### 3.1.2 CN Parameter Optimization
- [ ] **Graph size optimization**
  - CN detection: k = 20 neighbors
  - SC detection: k = 40 neighbors
  - Biological relevance validation

### 3.2 CN Composition and Characterization

#### 3.2.1 CN Composition Analysis
- [ ] **Cell type fraction heatmaps**
  - Raw counts per CN
  - Proportional composition
  - Column-scaled enrichment analysis
- [ ] **Statistical enrichment testing**
  - CN-celltype association testing
  - Multiple testing correction
  - Effect size quantification

#### 3.2.2 CN Spatial Distribution
- [ ] **Regional CN analysis**
  - CN abundance across tissues (CT/IM/TAT)
  - Spatial CN mapping with `hatchingPlot()`
  - CN boundary visualization

### 3.3 Clinical CN Correlation

#### 3.3.1 Clinical Group Comparison and Survival Analysis
- [ ] **Recurrence status comparison**
  - CN fraction differences (R vs NR)
  - Statistical testing with multiple comparisons
  - Effect size and clinical significance
  - 
- [ ] **CN abundance and survival**
  - Kaplan-Meier analysis per CN type
  - Cox regression modeling
  - Multivariable survival models
  - 
- [ ] **Treatment stratification**
  - CN patterns in chemotherapy vs combination therapy
  - Treatment-specific CN biomarkers
  - Treatment-CN interaction effects

### 3.4 CN Communication Analysis

#### 3.4.1 Canonical Correlation Analysis (CCA)
- [ ] **CN-CN correlation networks**
  - Pairwise CN correlation matrices
  - Network graph construction
  - Community detection in CN networks
- [ ] **Clinical group comparison**
  - CCA network differences (R vs NR)
  - Network topology analysis
  - Hub CN identification

#### 3.4.2 CN Interaction Networks
- [ ] **Spatial CN co-occurrence**
  - CN proximity analysis
  - Co-occurrence statistical testing
  - Network visualization

### 3.5 CN Functional Analysis

#### 3.5.1 CN-Specific Expression Analysis
- [ ] **Marker expression within CNs**
  - Mean expression per CN per patient
  - CN-specific functional signatures
  - Cross-CN expression comparison
- [ ] **Functional pathway analysis**
  - CN-specific pathway enrichment
  - Immune exhaustion signatures in immune CNs
  - Metabolic signatures in tumor CNs

#### 3.5.2 CN Functional Validation
- [ ] **scRNA-seq integration**
  - CN functional annotation validation
  - Gene expression program mapping
  - Single-cell CN composition validation

---

## 4. Advanced Spatial Analysis

### 4.1 Spatial Context (SC) Analysis

#### 4.1.1 SC Detection
- [ ] **`detectSpatialContext()` implementation**
  - CN fraction aggregation (k = 40)
  - Threshold optimization (0.90 default)
  - SC filtering by frequency and size
- [ ] **SC characterization**
  - SC composition analysis
  - Interface CN identification
  - SC network construction

#### 4.1.2 SC Clinical Correlation
- [ ] **SC-outcome association**
  - SC abundance and recurrence
  - SC-specific survival analysis
  - Treatment-SC interaction effects
- [ ] **SC network analysis**
  - `plotSpatialContext()` network graphs
  - SC connectivity analysis
  - Hub SC identification

### 4.2 Supervised Patch Detection

#### 4.2.1 Tumor Patch Analysis
- [ ] **`patchDetection()` implementation**
  - Connected tumor cell identification
  - Patch size distribution analysis
  - Tumor patch heterogeneity quantification
- [ ] **Invasion patch detection**
  - EMT-positive cell patches
  - Tumor budding patch identification
  - Invasion front characterization

#### 4.2.2 Microenvironment Patch Analysis
- [ ] **Immune exclusion patches**
  - Immune-depleted region identification
  - Exclusion patch size and frequency
  - Clinical correlation with outcomes
- [ ] **Metabolic patches**
  - Hypoxic cell cluster detection
  - Glycolytic patch identification
  - Metabolic heterogeneity analysis

### 4.3 Distance-Based Spatial Metrics

#### 4.3.1 Distance Analysis Implementation
- [ ] **`minDistToCells()` analysis**
  - Distance to tumor patches
  - Distance to immune cells
  - Distance to vasculature
- [ ] **Spatial gradient analysis**
  - Marker expression gradients
  - Distance-dependent expression analysis
  - Boundary effect quantification

#### 4.3.2 Clinical Distance Metrics
- [ ] **Immune infiltration quantification**
  - T cell distance to tumor cells
  - Macrophage proximity analysis
  - Immune exclusion quantification
- [ ] **Invasion distance metrics**
  - Distance from tumor core to margin
  - EMT gradient analysis
  - Budding cell distance analysis

---

## 5. Multi-Scale Spatial Integration

### 5.1 Scale-Dependent Analysis
- [ ] **Cellular scale** (5-10 μm)
  - Direct cell-cell interactions
  - Immediate neighborhood effects
- [ ] **Microenvironmental scale** (20-50 μm)
  - CN detection and characterization
  - Local tissue organization
- [ ] **Architectural scale** (100-200 μm)
  - SC analysis and patch detection
  - Regional organization patterns
- [ ] **Tissue scale** (mm)
  - Cross-regional comparison (CT/IM/TAT)
  - Global tissue architecture

### 5.2 Cross-Scale Integration
- [ ] **Multi-scale signature development**
  - Scale-specific feature extraction
  - Hierarchical clustering across scales
  - Integrated prognostic signatures

---

## 6. Advanced Statistical and Spatial Methods

### 6.1 Spatial Statistics
- [ ] **Spatial autocorrelation analysis**
  - Moran's I for expression gradients

### 6.2 Robustness and Validation
- [ ] **CN stability analysis**
  - Bootstrap resampling
  - Parameter sensitivity analysis
  - Adjusted Rand Index validation
- [ ] **Permutation testing**
  - Null model generation
  - Statistical significance validation
  - Multiple testing correction

---

---

## 10. Reproducibility and Documentation

### 10.1 Analysis Documentation
- [ ] **Parameter documentation**
  - All analysis parameters recorded
  - Rationale for parameter choices
  - Sensitivity analysis results

### 10.2 Code Reproducibility  
- [ ] **Version control**
  - All analysis code versioned
  - Environment documentation
  - Reproducible workflows

### 10.3 Result Validation
- [ ] **Independent validation**
  - External cohort testing
  - Cross-institutional validation
  - Clinical trial integration

---

## Expected Deliverables

### Scientific Outputs
- [ ] Comprehensive spatial interaction atlas
- [ ] CN-based prognostic signatures  
- [ ] Treatment response prediction models
- [ ] "Seed and soil" quantitative framework

### Clinical Outputs
- [ ] H&E-based diagnostic tool
- [ ] Risk stratification algorithm
- [ ] Treatment selection biomarkers
- [ ] Clinical implementation guidelines

### Technical Contributions
- [ ] Multi-modal spatial analysis pipeline
- [ ] Cross-platform validation framework
- [ ] Clinical translation methodology
- [ ] Reproducible analysis workflows