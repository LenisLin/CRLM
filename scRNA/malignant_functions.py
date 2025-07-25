# Comprehensive Clinical Analysis of NMF Clusters
import os
import copy
import pickle
import warnings
from itertools import combinations
from collections import Counter, defaultdict

warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import networkx as nx

import scanpy as sc
import gseapy as gp
from gseapy import Msigdb

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

import scipy
import scipy.io
import scipy.sparse
from scipy import stats
from scipy.stats import hypergeom, zscore, pearsonr
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform, pdist

from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import NMF, PCA

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def get_metabolism_pathways():
    """Enhanced metabolism pathway collection"""
    
    metabolism_gene_sets = {}
    
    # Enhanced metabolism keywords
    metabolism_keywords = [
        'metabolism', 'metabolic', 'glycolysis', 'glycolytic',
        'oxidative phosphorylation', 'electron transport', 'atp synthesis',
        'fatty acid', 'lipid', 'cholesterol', 'steroid',
        'amino acid', 'protein', 'nucleotide', 'purine', 'pyrimidine',
        'tca cycle', 'citrate cycle', 'krebs cycle',
        'pentose phosphate', 'hexosamine', 'fructose',
        'gluconeogenesis', 'glucose', 'insulin signaling',
        'autophagy', 'lysosome', 'peroxisome',
        'mitochondrial', 'ribosome', 'translation'
    ]
    
    # 1. KEGG Enhanced
    # try:
    #     kegg = gp.get_library(name='KEGG_2021_Human', organism='Human')
    #     kegg_metabolism = {k: v for k, v in kegg.items() 
    #                       if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     metabolism_gene_sets.update({f"KEGG_{k}": v for k, v in kegg_metabolism.items()})
    # except:
    #     pass
    
    # # 2. Reactome Enhanced
    # try:
    #     reactome = gp.get_library(name='Reactome_2022', organism='Human')
    #     reactome_metabolism = {k: v for k, v in reactome.items() 
    #                           if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     metabolism_gene_sets.update({f"REACTOME_{k}": v for k, v in reactome_metabolism.items()})
    # except:
    #     pass

    # 3. Hallmark Metabolism Pathways
    try:
        msig = Msigdb() # hallmark gene sets
        # msig.list_dbver() # list msigdb version you wanna query
        # msig.list_category(dbver="2025.1.Hs") # list categories given dbver.
        category = 'h.all'  # Hallmark gene sets
        dbver = "2025.1.Hs"  # Specify the version
        print(f"Fetching Hallmark gene sets from {category} with version {dbver}...")

        hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        hallmark_metabolism = {k: v for k, v in hallmark.items() 
                              if any(keyword in k.lower() for keyword in metabolism_keywords)}
        
        metabolism_gene_sets.update({f"HALLMARK_{k}": v for k, v in hallmark_metabolism.items()})
        print(f"✅ Hallmark metabolism pathways: {len(hallmark_metabolism)}")
        
    except Exception as e:
        print(f"⚠️ Hallmark access failed: {e}")

    
    # # 4. GO Biological Process - Metabolism
    # try:
    #     go_bp = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
    #     go_metabolism = {k: v for k, v in go_bp.items() 
    #                     if any(keyword in k.lower() for keyword in metabolism_keywords)}
    #     # Limit GO terms to prevent overload
    #     go_metabolism_top = dict(list(go_metabolism.items())[:100])
    #     metabolism_gene_sets.update({f"GO_BP_{k}": v for k, v in go_metabolism_top.items()})
    # except:
    #     pass
    
    return metabolism_gene_sets

def get_oncogenic_pathways():
    """Get oncogenic and cancer hallmark pathways"""
    
    oncogenic_gene_sets = {}
    
    # Cancer hallmark keywords
    cancer_keywords = [
        'cancer', 'tumor', 'oncogene', 'tumor suppressor',
        'cell cycle', 'apoptosis', 'proliferation', 'growth',
        'dna repair', 'dna damage', 'genome instability',
        'angiogenesis', 'invasion', 'metastasis',
        'emt', 'epithelial mesenchymal', 'stemness',
        'p53', 'rb', 'myc', 'ras', 'pi3k', 'akt', 'mtor',
        'wnt', 'notch', 'hedgehog', 'tgf', 'nf-kb'
    ]
    
    # 1. Hallmark Gene Sets
    try:
        hallmark = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        oncogenic_gene_sets.update({f"HALLMARK_{k}": v for k, v in hallmark.items()})
    except:
        pass
    
    # # 2. Oncogenic Signatures
    # try:
    #     oncogenic = gp.get_library(name='MSigDB_Oncogenic_Signatures', organism='Human')
    #     oncogenic_gene_sets.update({f"ONCOGENIC_{k}": v for k, v in oncogenic.items()})
    # except:
    #     pass
    
    # # 3. KEGG Cancer Pathways
    # try:
    #     kegg = gp.get_library(name='KEGG_2021_Human', organism='Human')
    #     kegg_cancer = {k: v for k, v in kegg.items() 
    #                   if any(keyword in k.lower() for keyword in cancer_keywords)}
    #     oncogenic_gene_sets.update({f"KEGG_CANCER_{k}": v for k, v in kegg_cancer.items()})
    # except:
    #     pass
    
    return oncogenic_gene_sets

def get_immune_microenvironment_pathways():
    """Get immune and microenvironment pathways"""
    
    immune_gene_sets = {}
    
    immune_keywords = [
        'immune', 'immunity', 'immunological',
        'interferon', 'interleukin', 'cytokine', 'chemokine',
        'inflammatory', 'inflammation', 'response',
        't cell', 'b cell', 'nk cell', 'macrophage', 'dendritic',
        'antigen presentation', 'mhc', 'hla',
        'complement', 'toll like', 'innate immunity',
        'adaptive immunity', 'antibody', 'immunoglobulin'
    ]
    
    try:
        # Reactome immune pathways
        reactome = gp.get_library(name='Reactome_2022', organism='Human')
        reactome_immune = {k: v for k, v in reactome.items() 
                          if any(keyword in k.lower() for keyword in immune_keywords)}
        immune_gene_sets.update({f"REACTOME_IMMUNE_{k}": v for k, v in reactome_immune.items()})
        
        # GO immune process
        go_bp = gp.get_library(name='GO_Biological_Process_2023', organism='Human')
        go_immune = {k: v for k, v in go_bp.items() 
                    if any(keyword in k.lower() for keyword in immune_keywords)}
        go_immune_top = dict(list(go_immune.items())[:50])
        immune_gene_sets.update({f"GO_IMMUNE_{k}": v for k, v in go_immune_top.items()})
        
    except:
        pass
    
    return immune_gene_sets

def get_stress_resistance_pathways():
    """Get stress response and treatment resistance pathways"""
    
    stress_gene_sets = {}
    
    stress_keywords = [
        'stress', 'response', 'resistance', 'drug resistance',
        'hypoxia', 'oxidative stress', 'er stress', 'unfolded protein',
        'heat shock', 'chaperone', 'proteostasis',
        'autophagy', 'apoptosis', 'necroptosis', 'ferroptosis',
        'dna damage response', 'checkpoint', 'repair',
        'chemotherapy', 'radiation', 'therapy resistance'
    ]
    
    try:
        # Hallmark stress pathways
        hallmark = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        hallmark_stress = {k: v for k, v in hallmark.items() 
                          if any(keyword in k.lower() for keyword in stress_keywords)}
        stress_gene_sets.update({f"HALLMARK_STRESS_{k}": v for k, v in hallmark_stress.items()})
        
        # # Reactome stress pathways
        # reactome = gp.get_library(name='Reactome_2022', organism='Human')
        # reactome_stress = {k: v for k, v in reactome.items() 
        #                   if any(keyword in k.lower() for keyword in stress_keywords)}
        # stress_gene_sets.update({f"REACTOME_STRESS_{k}": v for k, v in reactome_stress.items()})
        
    except:
        pass
    
    return stress_gene_sets

def get_highly_variable_genes(adata, candidate_genes = None, n_top_genes=3000):
    """Get highly variable genes from candidates"""
    
    if candidate_genes:
        # First filter to only include candidate genes that exist in the data
        available_candidates = [gene for gene in candidate_genes if gene in adata.var_names]
        
        if len(available_candidates) == 0:
            print("Warning: No candidate genes found in the data")
            return []
        
        # Subset the data to only include candidate genes
        adata_subset = adata[:, available_candidates].copy()

        # Adjust n_top_genes if we have fewer candidates than requested
        n_genes_to_select = min(n_top_genes, len(available_candidates))
        
    else:
        adata_subset = adata.copy()
        n_genes_to_select = min(n_top_genes, adata_subset.shape[1])

    # Perform HVG selection on the filtered subset
    sc.pp.highly_variable_genes(adata_subset, n_top_genes=n_genes_to_select)
    
    # Get the highly variable genes from the filtered set
    hvg_genes = adata_subset.var_names[adata_subset.var.highly_variable].tolist()
        
    return hvg_genes

def get_multi_pathway_genes(gene_sets, candidate_genes, min_pathways=2):
    """Get genes that appear in multiple pathways"""
    
    gene_pathway_count = defaultdict(int)
    
    for pathway, genes in gene_sets.items():
        for gene in genes:
            if gene in candidate_genes:
                gene_pathway_count[gene] += 1
    
    multi_pathway_genes = [gene for gene, count in gene_pathway_count.items() 
                          if count >= min_pathways]
    
    return multi_pathway_genes

def get_cancer_associated_genes(candidate_genes):
    """Get known cancer-associated genes"""
    
    # Common cancer genes (you can expand this list)
    cancer_genes = {
        'TP53', 'KRAS', 'PIK3CA', 'APC', 'PTEN', 'EGFR', 'MYC', 'RB1',
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'MLH1', 'MSH2', 'MSH6', 'PMS2',
        'VHL', 'NF1', 'NF2', 'CDKN2A', 'CDK4', 'MDM2', 'ERBB2', 'MET',
        'ALK', 'ROS1', 'RET', 'NTRK1', 'BRAF', 'NRAS', 'HRAS', 'KIT',
        'PDGFRA', 'FLT3', 'IDH1', 'IDH2', 'TERT', 'ARID1A', 'CTNNB1'
    }
    
    cancer_candidates = [gene for gene in candidate_genes if gene in cancer_genes]
    
    return cancer_candidates

def get_imc_aligned_hallmark_pathways():
    """
    Get HALLMARK pathways specifically aligned with IMC panel markers
    Focused on pathways detectable by: GLUT1, HK2, FASN, PRPS1, CA-IX, VEGF, VIM, EPCAM, Ki67
    """
    
    imc_aligned_gene_sets = {}
    
    try:
        msig = Msigdb()
        
        print("🎯 Fetching IMC-aligned HALLMARK pathways...")
        
        # Get all HALLMARK gene sets
        hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        
        # Define IMC-aligned HALLMARK pathways
        imc_relevant_hallmarks = {
            # Metabolism pathways (align with GLUT1, HK2, FASN)
            'HALLMARK_GLYCOLYSIS': 'Genes involved in glycolysis (→ GLUT1, HK2)',
            'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'Genes encoding oxidative phosphorylation machinery',
            'HALLMARK_FATTY_ACID_METABOLISM': 'Genes regulating fatty acid metabolism (→ FASN)',
            
            # Proliferation pathways (align with Ki67)
            'HALLMARK_E2F_TARGETS': 'Genes encoding E2F transcription factor targets (→ Ki67)',
            'HALLMARK_G2M_CHECKPOINT': 'Genes involved in G2/M checkpoint (→ Ki67)',
            'HALLMARK_MYC_TARGETS_V1': 'Genes encoding MYC targets (proliferation)',
            'HALLMARK_MYC_TARGETS_V2': 'Genes encoding MYC targets (metabolic)',
            
            # Stress/Hypoxia pathways (align with CA-IX, VEGF)
            'HALLMARK_HYPOXIA': 'Genes up-regulated in hypoxia (→ CA-IX, VEGF)',
            'HALLMARK_ANGIOGENESIS': 'Genes involved in angiogenesis (→ VEGF)',
            
            # EMT/Invasion pathways (align with Vimentin, EpCAM)
            'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': 'Genes defining EMT (→ VIM, EPCAM)',
            
            # Growth/Anabolic pathways
            'HALLMARK_MTORC1_SIGNALING': 'Genes up-regulated by mTORC1 activation',
            
            # DNA damage/Cell cycle control
            'HALLMARK_P53_PATHWAY': 'Genes involved in p53 pathway and DNA damage response',
            'HALLMARK_DNA_REPAIR': 'Genes involved in DNA repair',
            
            # Apoptosis (complement to proliferation)
            'HALLMARK_APOPTOSIS': 'Genes mediating programmed cell death'
        }
        
        # Extract only the IMC-relevant pathways
        for pathway_name, description in imc_relevant_hallmarks.items():
            if pathway_name in hallmark:
                imc_aligned_gene_sets[pathway_name] = hallmark[pathway_name]
                print(f"✅ {pathway_name}: {len(hallmark[pathway_name])} genes")
            else:
                print(f"⚠️ {pathway_name}: Not found in HALLMARK database")
        
        print(f"\n🎯 Total IMC-aligned HALLMARK pathways: {len(imc_aligned_gene_sets)}")
        
        # Print pathway-IMC marker mapping for verification
        print("\n📊 Pathway → IMC Marker Alignment:")
        pathway_imc_mapping = {
            'HALLMARK_GLYCOLYSIS': ['GLUT1', 'HK2'],
            'HALLMARK_FATTY_ACID_METABOLISM': ['FASN'],
            'HALLMARK_HYPOXIA': ['CA-IX', 'VEGF'],
            'HALLMARK_ANGIOGENESIS': ['VEGF'],
            'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': ['Vimentin', 'EpCAM(inverse)'],
            'HALLMARK_E2F_TARGETS': ['Ki67'],
            'HALLMARK_G2M_CHECKPOINT': ['Ki67'],
            'HALLMARK_MYC_TARGETS_V1': ['Ki67', 'metabolic markers'],
            'HALLMARK_MYC_TARGETS_V2': ['FASN', 'metabolic markers']
        }
        
        for pathway, imc_markers in pathway_imc_mapping.items():
            if pathway in imc_aligned_gene_sets:
                print(f"   {pathway} → {', '.join(imc_markers)}")
        
    except Exception as e:
        print(f"❌ Error fetching HALLMARK pathways: {e}")
        
        # Fallback: Manual definition of key pathways (if MSigDB fails)
        print("🔄 Using fallback manual pathway definitions...")
    
    return imc_aligned_gene_sets

def get_all_hallmark_pathways():
    """
    Get HALLMARK pathways specifically aligned with IMC panel markers
    Focused on pathways detectable by: GLUT1, HK2, FASN, PRPS1, CA-IX, VEGF, VIM, EPCAM, Ki67
    """
        
    msig = Msigdb()
    
    print("🎯 Fetching IMC-aligned HALLMARK pathways...")
    
    # Get all HALLMARK gene sets
    hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
    
    return hallmark

def validate_imc_pathway_coverage(pathways_dict, malignant_cells):
    """
    Validate how well the pathways cover IMC panel genes
    """
    
    print("\n🔍 Validating IMC panel coverage in pathways...")
    
    # IMC panel gene symbols
    imc_panel_genes = {
        'SLC2A1': 'GLUT1',      # Glycolysis
        'HK2': 'HK2',           # Glycolysis  
        'FASN': 'FASN',         # Fatty acid metabolism
        'PRPS1': 'PRPS1',       # Nucleotide metabolism (may not be in HALLMARK)
        'CA9': 'CA-IX',         # Hypoxia
        'VEGFA': 'VEGF',        # Angiogenesis/Hypoxia
        'VIM': 'Vimentin',      # EMT
        'EPCAM': 'EpCAM',       # Epithelial (inverse of EMT)
        'MKI67': 'Ki67'         # Proliferation
    }
    
    # Check coverage
    coverage_report = {}
    
    for gene_symbol, protein_name in imc_panel_genes.items():
        found_in_pathways = []
        
        for pathway_name, pathway_genes in pathways_dict.items():
            if gene_symbol in pathway_genes:
                found_in_pathways.append(pathway_name)
        
        # Check if gene exists in dataset
        in_dataset = gene_symbol in malignant_cells.var_names
        
        coverage_report[protein_name] = {
            'gene_symbol': gene_symbol,
            'in_dataset': in_dataset,
            'in_pathways': found_in_pathways,
            'n_pathways': len(found_in_pathways)
        }
    
    # Print coverage report
    print("\n📊 IMC Panel Gene Coverage:")
    print("=" * 60)
    
    for protein, info in coverage_report.items():
        status = "✅" if info['in_dataset'] else "❌"
        pathway_count = info['n_pathways']
        
        print(f"{status} {protein} ({info['gene_symbol']}): "
              f"{pathway_count} pathways")
        
        if info['in_pathways']:
            for pathway in info['in_pathways']:
                print(f"    - {pathway}")
    
    return coverage_report

def get_overlap_count(list1, list2):
    return len(set(list1) & set(list2))

def calculate_gene_overlap(list1, list2):
    intersection = len(set(list1) & set(list2))
    union = len(set(list1) | set(list2))
    return intersection / union if union > 0 else 0

# def run_robust_nmf_analysis(malignant_cells, available_metabolism_genes, save_path):
#     """
#     Implement robust NMF strategy following CSCC paper methodology
#     """
    
#     print("🎯 Starting robust NMF analysis...")
    
#     # Parameters following the paper
#     K_RANGE = range(4, 9)  # Try: range(3, 8) or range(4, 10)
#     TOP_GENES_PER_PROGRAM = 50  # Top 25 genes per program
#     MIN_CELLS_PER_PATIENT = 150  # Minimum 150 cells per patient
#     MAX_CLUSTER_NUM = 8

#     # Robustness criteria
#     SAME_PATIENT_OVERLAP_THRESHOLD = 0.7  # 50% overlap
#     CROSS_PATIENT_OVERLAP_THRESHOLD = 0.2  # 20% overlap
#     MAX_INTRA_PATIENT_OVERLAP = 0.2  # <20% overlap within same patient
    
#     # Group cells by patient/study
#     patient_groups = malignant_cells.obs['sample_id'].unique()
#     print(f"Found {len(patient_groups)} unique samples")
    
#     # Filter patients with sufficient cells
#     valid_patients = []
#     for patient in patient_groups:
#         patient_cells = malignant_cells[malignant_cells.obs['sample_id'] == patient]
#         if patient_cells.n_obs >= MIN_CELLS_PER_PATIENT:
#             valid_patients.append(patient)
    
#     print(f"{len(valid_patients)} samples with ≥{MIN_CELLS_PER_PATIENT} cells")
    
#     # Step 1: Run NMF for each patient and each K
#     print("\n📊 Step 1: Running NMF for each patient...")
#     all_programs = []
#     patient_programs = {}
    
#     for patient_id in valid_patients:
#         print(f"  Processing sample {patient_id}...")
        
#         # Subset patient cells
#         patient_cells = malignant_cells[malignant_cells.obs['sample_id'] == patient_id]
#         patient_metabolism = patient_cells[:, available_metabolism_genes].copy()
        
#         # Preprocessing
#         # sc.pp.normalize_total(patient_metabolism, target_sum=1e4)
#         if scipy.sparse.issparse(patient_metabolism.X):
#             X_patient = patient_metabolism.X.toarray()
#         else:
#             X_patient = patient_metabolism.X.copy()
        
#         X_patient = np.maximum(X_patient, 0)  # Ensure non-negative
        
#         patient_programs[patient_id] = {}
        
#         # Run NMF for different K values
#         for k in K_RANGE:
#             print(f"    K={k}: {patient_metabolism.n_obs} cells, {patient_metabolism.n_vars} genes")
            
#             # Run NMF
#             nmf_model = NMF(n_components=k, random_state=619, max_iter=2000)
#             W = nmf_model.fit_transform(X_patient)
#             H = nmf_model.components_
            
#             # Extract top genes for each program
#             k_programs = []
#             for factor_idx in range(k):
#                 # Get top 50 genes for this factor
#                 factor_scores = H[factor_idx, :]
#                 top_gene_indices = np.argsort(factor_scores)[-TOP_GENES_PER_PROGRAM:][::-1]
#                 top_genes = [available_metabolism_genes[i] for i in top_gene_indices]
                
#                 program = {
#                     'patient': patient_id,
#                     'k_value': k,
#                     'factor_idx': factor_idx,
#                     'genes': top_genes,
#                     'gene_scores': factor_scores[top_gene_indices],
#                     'program_id': f"{patient_id}_K{k}_F{factor_idx}"
#                 }
                
#                 k_programs.append(program)
#                 all_programs.append(program)
            
#             patient_programs[patient_id][k] = k_programs
    
#     print(f"✅ Generated {len(all_programs)} total programs from {len(valid_patients)} patients")
    
#     # Step 2: Apply robustness criteria
#     print("\n🔍 Step 2: Applying robustness criteria...")
    
#     # Criterion 1: 70% overlap with different K in same patient
#     print("  Applying criterion 1: Cross-K stability within patients...")
#     stable_programs = []

#     for program in all_programs:
#         patient_id = program['patient']
#         k_value = program['k_value']
#         genes = program['genes']
        
#         # Check overlap with other K values for same patient
#         cross_k_overlaps = []
#         for other_k in K_RANGE:
#             if other_k != k_value and other_k in patient_programs[patient_id]:
#                 for other_program in patient_programs[patient_id][other_k]:
#                     overlap = get_overlap_count(genes, other_program['genes'])
#                     cross_k_overlaps.append(overlap)
        
#         # Check if any overlap meets the 70% threshold (35/50 genes)
#         same_paient_threshold = int(TOP_GENES_PER_PROGRAM * SAME_PATIENT_OVERLAP_THRESHOLD)

#         if any(overlap >= same_paient_threshold for overlap in cross_k_overlaps):
#             stable_programs.append(program)
    
#     print(f"    Programs passing criterion 1: {len(stable_programs)}")
    
#     # Criterion 2: 20% overlap with programs in other patients
#     print("  Applying criterion 2: Cross-patient similarity...")
#     cross_patient_programs = []

#     cross_patient_threshold = int(TOP_GENES_PER_PROGRAM * CROSS_PATIENT_OVERLAP_THRESHOLD)
    
#     for program in stable_programs:
#         patient_id = program['patient']
#         genes = program['genes']
        
#         # Check overlap with programs from other patients
#         cross_patient_overlaps = []
#         for other_program in stable_programs:
#             if other_program['patient'] != patient_id:
#                 overlap = get_overlap_count(genes, other_program['genes'])
#                 cross_patient_overlaps.append(overlap)
        
#         # Check if any overlap meets the 20% threshold (10/50 genes)
#         if any(overlap >= cross_patient_threshold for overlap in cross_patient_overlaps):
#             cross_patient_programs.append(program)
    
#     print(f"    Programs passing criterion 2: {len(cross_patient_programs)}")
    
#     # Criterion 3: <20% overlap with other programs in same patient + selection
#     print("  Applying criterion 3: Intra-patient redundancy removal...")
    
#     # Rank programs by cross-patient similarity
#     program_similarities = []
#     for program in cross_patient_programs:
#         patient_id = program['patient']
#         genes = program['genes']
        
#         # Calculate average similarity to programs in other patients
#         similarities = []
#         for other_program in cross_patient_programs:
#             if other_program['patient'] != patient_id:
#                 overlap = get_overlap_count(genes, other_program['genes'])
#                 similarities.append(overlap)
        
#         avg_similarity = np.mean(similarities) if similarities else 0
#         program_similarities.append((program, avg_similarity))
    
#     # Sort by decreasing cross-patient similarity
#     program_similarities.sort(key=lambda x: x[1], reverse=True)
    
#     # Select programs, avoiding redundancy within patients
#     robust_programs = []
#     selected_patient_programs = {patient: [] for patient in valid_patients}
    
#     high_overlap_threshold = int(TOP_GENES_PER_PROGRAM * MAX_INTRA_PATIENT_OVERLAP)

#     for program, similarity in program_similarities:
#         patient_id = program['patient']
#         genes = program['genes']
        
#         # Check overlap with already selected programs from same patient
#         intra_patient_overlaps = []
#         for selected_program in selected_patient_programs[patient_id]:
#             overlap = get_overlap_count(genes, selected_program['genes'])
#             intra_patient_overlaps.append(overlap)
        
#         # Select if no high overlap with selected programs from same patient
#         if not any(overlap >= high_overlap_threshold for overlap in intra_patient_overlaps):
#             robust_programs.append(program)
#             selected_patient_programs[patient_id].append(program)
    
#     print(f"    Final robust programs: {len(robust_programs)}")
    
#     # Step 3: Hierarchical clustering to identify Meta-Programs (MPs)
#     print("\n🌳 Step 3: Hierarchical clustering to identify Meta-Programs...")
    
#     if len(robust_programs) < 2:
#         print("  ⚠️ Too few robust programs for clustering")
#         return None
    
#     # Create similarity matrix using Jaccard similarity
#     n_programs = len(robust_programs)
#     similarity_matrix = np.zeros((n_programs, n_programs))
    
#     for i in range(n_programs):
#         for j in range(n_programs):
#             overlap = calculate_gene_overlap(robust_programs[i]['genes'], robust_programs[j]['genes'])
#             similarity_matrix[i, j] = overlap
    
#     # Convert to distance matrix
#     distance_matrix = 1 - similarity_matrix
    
#     # Hierarchical clustering
#     condensed_distances = squareform(distance_matrix, checks=False)
#     linkage_matrix = linkage(condensed_distances, method='average')
    
#     # Determine optimal number of clusters (Meta-Programs)
#     # Try different numbers of clusters and evaluate
#     max_clusters = min(MAX_CLUSTER_NUM, len(robust_programs))
#     cluster_solutions = {}
    
#     for n_clusters in range(2, max_clusters + 1):
#         clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        
#         # Calculate within-cluster similarity
#         within_cluster_similarities = []
#         for cluster_id in range(1, n_clusters + 1):
#             cluster_programs = [i for i, c in enumerate(clusters) if c == cluster_id]
#             if len(cluster_programs) > 1:
#                 cluster_similarities = []
#                 for i, j in combinations(cluster_programs, 2):
#                     cluster_similarities.append(similarity_matrix[i, j])
#                 within_cluster_similarities.extend(cluster_similarities)
        
#         avg_within_similarity = np.mean(within_cluster_similarities) if within_cluster_similarities else 0
#         cluster_solutions[n_clusters] = {
#             'clusters': clusters,
#             'within_similarity': avg_within_similarity,
#             'n_programs_per_cluster': [np.sum(clusters == i) for i in range(1, n_clusters + 1)]
#         }
    
#     # Select optimal clustering (highest within-cluster similarity)
#     optimal_n_clusters = max(cluster_solutions.keys(), 
#                            key=lambda k: cluster_solutions[k]['within_similarity'])
#     optimal_clusters = cluster_solutions[optimal_n_clusters]['clusters']
    
#     print(f"    Optimal number of Meta-Programs: {optimal_n_clusters}")
#     print(f"    Programs per MP: {cluster_solutions[optimal_n_clusters]['n_programs_per_cluster']}")
    
#     # Step 4: Generate Meta-Program signatures
#     print("\n🧬 Step 4: Generating Meta-Program signatures...")
    
#     meta_programs = {}
    
#     for mp_id in range(1, optimal_n_clusters + 1):
#         mp_program_indices = [i for i, c in enumerate(optimal_clusters) if c == mp_id]
#         mp_programs = [robust_programs[i] for i in mp_program_indices]
        
#         print(f"  Meta-Program {mp_id}: {len(mp_programs)} constituent programs")
        
#         # Collect all genes from constituent programs
#         all_mp_genes = []
#         for program in mp_programs:
#             all_mp_genes.extend(program['genes'])
        
#         # Count gene occurrences
#         gene_counts = Counter(all_mp_genes)
        
#         # Select genes occurring in ≥25% of constituent programs
#         min_occurrence = max(1, len(mp_programs) * 0.25)
#         candidate_genes = [gene for gene, count in gene_counts.items() if count >= min_occurrence]
        
#         # Calculate gene scores (average across programs where gene appears)
#         gene_scores = {}
#         for gene in candidate_genes:
#             scores = []
#             for program in mp_programs:
#                 if gene in program['genes']:
#                     gene_idx = program['genes'].index(gene)
#                     scores.append(program['gene_scores'][gene_idx])
#             gene_scores[gene] = np.mean(scores)
        
#         # Sort by score and take top genes
#         sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
        
#         # Limit to reasonable signature size
#         max_signature_size = min(50, len(sorted_genes))
#         signature_genes = [gene for gene, score in sorted_genes[:max_signature_size]]
#         signature_scores = [score for gene, score in sorted_genes[:max_signature_size]]
        
#         meta_programs[mp_id] = {
#             'mp_id': mp_id,
#             'constituent_programs': mp_programs,
#             'signature_genes': signature_genes,
#             'signature_scores': signature_scores,
#             'n_programs': len(mp_programs),
#             'patients_represented': list(set(p['patient'] for p in mp_programs))
#         }
        
#         print(f"    MP{mp_id} signature: {len(signature_genes)} genes, "
#               f"{len(meta_programs[mp_id]['patients_represented'])} samples")
    
#     # Step 5: Calculate Meta-Program scores for all cells
#     print("\n📊 Step 5: Calculating Meta-Program scores for all cells...")
    
#     # Calculate MP scores using scanpy
#     for mp_id, mp_data in meta_programs.items():
#         signature_genes = mp_data['signature_genes']
#         available_signature_genes = [g for g in signature_genes if g in malignant_cells.var_names]
        
#         if len(available_signature_genes) >= 5:  # Minimum genes for scoring
#             sc.tl.score_genes(malignant_cells, available_signature_genes, 
#                              score_name=f'MP{mp_id}_score', use_raw=False)
#             print(f"    MP{mp_id}: {len(available_signature_genes)} genes used for scoring")
    
#     # Step 6: Assign cells to Meta-Programs
#     print("\n🎯 Step 6: Assigning cells to Meta-Programs...")
    
#     # Get MP score columns
#     mp_score_cols = [f'MP{mp_id}_score' for mp_id in meta_programs.keys() 
#                      if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
#     if len(mp_score_cols) >= 2:
#         # Normalize scores by mean subtraction
#         for col in mp_score_cols:
#             malignant_cells.obs[f'{col}_normalized'] = (malignant_cells.obs[col] - 
#                                                        malignant_cells.obs[col].mean())
        
#         # Assign cells to MPs
#         normalized_cols = [f'{col}_normalized' for col in mp_score_cols]
#         mp_assignments = []
        
#         for idx in range(malignant_cells.n_obs):
#             scores = [malignant_cells.obs.iloc[idx][col] for col in normalized_cols]
            
#             if len(scores) >= 2:
#                 max_score_idx = np.argmax(scores)
#                 max_score = scores[max_score_idx]
                
#                 # Get second highest score
#                 scores_sorted = sorted(scores, reverse=True)
#                 second_highest = scores_sorted[1] if len(scores_sorted) > 1 else 0
                
#                 # Assign if highest score > 85% threshold of second highest
#                 if second_highest == 0 or max_score > 0.85 * abs(second_highest):
#                     mp_assignments.append(f'MP{max_score_idx + 1}')
#                 else:
#                     mp_assignments.append('Unresolved')
#             else:
#                 mp_assignments.append('Unresolved')
        
#         malignant_cells.obs['MP_assignment'] = mp_assignments
#         malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
        
#         # Print assignment statistics
#         assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
#         print("    MP assignment distribution:")
#         for mp, count in assignment_counts.items():
#             percentage = count / len(malignant_cells) * 100
#             print(f"      {mp}: {count} cells ({percentage:.1f}%)")
    
#     # Save results
#     print("\n💾 Saving robust NMF results...")
    
#     robust_nmf_results = {
#         'all_programs': all_programs,
#         'robust_programs': robust_programs,
#         'meta_programs': meta_programs,
#         'clustering_info': {
#             'optimal_n_clusters': optimal_n_clusters,
#             'similarity_matrix': similarity_matrix,
#             'linkage_matrix': linkage_matrix,
#             'cluster_assignments': optimal_clusters
#         },
#         'parameters': {
#             'K_range': list(K_RANGE),
#             'top_genes_per_program': TOP_GENES_PER_PROGRAM,
#             'min_cells_per_patient': MIN_CELLS_PER_PATIENT,
#             'same_patient_overlap_threshold': SAME_PATIENT_OVERLAP_THRESHOLD,
#             'cross_patient_overlap_threshold': CROSS_PATIENT_OVERLAP_THRESHOLD,
#             'max_intra_patient_overlap': MAX_INTRA_PATIENT_OVERLAP
#         }
#     }
    
#     with open(os.path.join(save_path, "robust_metabolism_nmf_results.pkl"), 'wb') as f:
#         pickle.dump(robust_nmf_results, f)
    
#     print(f"✅ Robust NMF analysis completed!")
#     print(f"📊 Summary:")
#     print(f"   - {len(all_programs)} total programs generated")
#     print(f"   - {len(robust_programs)} robust programs identified")  
#     print(f"   - {optimal_n_clusters} Meta-Programs discovered")
#     print(f"   - {len(valid_patients)} patients analyzed")
    
#     return robust_nmf_results

def run_robust_nmf_analysis_aligned(malignant_cells, available_metabolism_genes, save_path):
    """
    Implement robust NMF strategy aligned with the CSCC paper methodology.
    """
    
    print("🎯 Starting robust NMF analysis (Aligned Version)...")
    
    # Parameters following the paper
    K_RANGE = range(4, 10)
    TOP_GENES_PER_PROGRAM = 50
    MIN_CELLS_PER_PATIENT = 150
    N_META_PROGRAMS = 4

    # Robustness criteria
    SAME_PATIENT_OVERLAP_THRESHOLD = 0.7
    CROSS_PATIENT_OVERLAP_THRESHOLD = 0.2
    MAX_INTRA_PATIENT_OVERLAP = 0.2
    
    # Group cells by patient/study
    patient_groups = malignant_cells.obs['sample_id'].unique()
    print(f"Found {len(patient_groups)} unique samples")
    
    # Filter patients with sufficient cells
    valid_patients = [
        p for p in patient_groups 
        if malignant_cells[malignant_cells.obs['sample_id'] == p].n_obs >= MIN_CELLS_PER_PATIENT
    ]
    
    print(f"{len(valid_patients)} samples with ≥{MIN_CELLS_PER_PATIENT} cells")
    
    # Step 1: Run NMF for each patient and each K
    print("\n📊 Step 1: Running NMF for each patient...")
    all_programs = []
    patient_programs = {}
    
    for patient_id in valid_patients:
        print(f"  Processing sample {patient_id}...")
        
        patient_cells = malignant_cells[malignant_cells.obs['sample_id'] == patient_id]
        patient_metabolism = patient_cells[:, available_metabolism_genes].copy()
        
        if scipy.sparse.issparse(patient_metabolism.X):
            X_patient = patient_metabolism.X.toarray()
        else:
            X_patient = patient_metabolism.X.copy()
        
        X_patient = np.maximum(X_patient, 0)
        
        patient_programs[patient_id] = {}
        
        for k in K_RANGE:
            print(f"    K={k}: {patient_metabolism.n_obs} cells, {patient_metabolism.n_vars} genes")
            nmf_model = NMF(n_components=k, random_state=619, max_iter=2000, init='nndsvda')
            W = nmf_model.fit_transform(X_patient)
            H = nmf_model.components_
            
            k_programs = []
            for factor_idx in range(k):
                factor_scores = H[factor_idx, :]
                top_gene_indices = np.argsort(factor_scores)[-TOP_GENES_PER_PROGRAM:][::-1]
                top_genes = [available_metabolism_genes[i] for i in top_gene_indices]
                
                program = {
                    'patient': patient_id, 'k_value': k, 'factor_idx': factor_idx,
                    'genes': top_genes, 'gene_scores': factor_scores[top_gene_indices],
                    'program_id': f"{patient_id}_K{k}_F{factor_idx}"
                }
                k_programs.append(program)
                all_programs.append(program)
            
            patient_programs[patient_id][k] = k_programs
    
    print(f"✅ Generated {len(all_programs)} total programs from {len(valid_patients)} patients")
    
    # Step 2: Apply robustness criteria
    print("\n🔍 Step 2: Applying robustness criteria...")
    
    # Criterion 1: Cross-K stability within patients
    stable_programs = []
    same_patient_threshold = int(TOP_GENES_PER_PROGRAM * SAME_PATIENT_OVERLAP_THRESHOLD)
    for program in all_programs:
        patient_id, k_value, genes = program['patient'], program['k_value'], program['genes']
        has_stable_overlap = False
        for other_k in K_RANGE:
            if other_k != k_value:
                for other_program in patient_programs[patient_id].get(other_k, []):
                    if get_overlap_count(genes, other_program['genes']) >= same_patient_threshold:
                        has_stable_overlap = True
                        break
            if has_stable_overlap:
                break
        if has_stable_overlap:
            stable_programs.append(program)
    print(f"  Programs passing criterion 1: {len(stable_programs)}")

    # Criterion 2: Cross-patient similarity
    cross_patient_programs = []
    cross_patient_threshold = int(TOP_GENES_PER_PROGRAM * CROSS_PATIENT_OVERLAP_THRESHOLD)
    for program in stable_programs:
        has_cross_overlap = False
        for other_program in stable_programs:
            if other_program['patient'] != program['patient']:
                if get_overlap_count(program['genes'], other_program['genes']) >= cross_patient_threshold:
                    has_cross_overlap = True
                    break
        if has_cross_overlap:
            cross_patient_programs.append(program)
    print(f"  Programs passing criterion 2: {len(cross_patient_programs)}")

    # Criterion 3: Intra-patient redundancy removal
    program_similarities = []
    for program in cross_patient_programs:
        similarities = [
            get_overlap_count(program['genes'], other['genes'])
            for other in cross_patient_programs if other['patient'] != program['patient']
        ]
        avg_similarity = np.mean(similarities) if similarities else 0
        program_similarities.append((program, avg_similarity))
    
    program_similarities.sort(key=lambda x: x[1], reverse=True)
    
    robust_programs = []
    selected_in_patient = {p: [] for p in valid_patients}
    high_overlap_threshold = int(TOP_GENES_PER_PROGRAM * MAX_INTRA_PATIENT_OVERLAP)
    
    for program, _ in program_similarities:
        patient_id = program['patient']
        is_redundant = False
        for selected in selected_in_patient[patient_id]:
            if get_overlap_count(program['genes'], selected['genes']) >= high_overlap_threshold:
                is_redundant = True
                break
        if not is_redundant:
            robust_programs.append(program)
            selected_in_patient[patient_id].append(program)
    print(f"  Final robust programs: {len(robust_programs)}")
    
    # Step 3: Hierarchical clustering to identify Meta-Programs (MPs)
    print("\n🌳 Step 3: Hierarchical clustering to identify Meta-Programs...")
    
    if len(robust_programs) < N_META_PROGRAMS:
        print(f"  ⚠️ Too few robust programs ({len(robust_programs)}) for clustering into {N_META_PROGRAMS} MPs")
        return None
        
    n_programs = len(robust_programs)
    similarity_matrix = np.array([[calculate_gene_overlap(p1['genes'], p2['genes']) for p2 in robust_programs] for p1 in robust_programs])
    distance_matrix = 1 - similarity_matrix
    
    condensed_distances = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_distances, method='average') 
    
    optimal_n_clusters = N_META_PROGRAMS
    clusters = fcluster(linkage_matrix, optimal_n_clusters, criterion='maxclust')
    print(f"  Clustered into {optimal_n_clusters} Meta-Programs (fixed)")
    
    # Step 4: Generate Meta-Program signatures
    print("\n🧬 Step 4: Generating Meta-Program signatures...")
    
    meta_programs = {}
    
    full_metabolism_adata = malignant_cells[:, available_metabolism_genes]
    if scipy.sparse.issparse(full_metabolism_adata.X):
        full_expr_matrix = full_metabolism_adata.X.toarray()
    else:
        full_expr_matrix = full_metabolism_adata.X
    
    expr_df = pd.DataFrame(full_expr_matrix, columns=available_metabolism_genes)

    for mp_id in range(1, optimal_n_clusters + 1):
        mp_program_indices = [i for i, c in enumerate(clusters) if c == mp_id]
        if not mp_program_indices:
            continue
            
        mp_programs = [robust_programs[i] for i in mp_program_indices]
        print(f"  Meta-Program {mp_id}: {len(mp_programs)} constituent programs")
        
        all_mp_genes = [gene for prog in mp_programs for gene in prog['genes']]
        gene_counts = Counter(all_mp_genes)
        min_occurrence = max(1, len(mp_programs) * 0.25)
        candidate_genes = [gene for gene, count in gene_counts.items() if count >= min_occurrence]
        
        if len(candidate_genes) < 2:
            print(f"    ⚠️ Not enough candidate genes for MP {mp_id}, skipping.")
            continue

        temp_score_name = f'temp_MP{mp_id}_score'
        sc.tl.score_genes(malignant_cells, candidate_genes, score_name=temp_score_name, use_raw=False)
        mp_signature_scores = malignant_cells.obs[temp_score_name].values
        
        correlations = expr_df.corrwith(pd.Series(mp_signature_scores, index=expr_df.index), method='pearson')
        top_30_correlated = correlations.nlargest(30)
        
        signature_genes = top_30_correlated.index.tolist()
        signature_scores = top_30_correlated.values.tolist()

        del malignant_cells.obs[temp_score_name]
        
        meta_programs[mp_id] = {
            'mp_id': mp_id,
            'constituent_programs': mp_programs,
            'signature_genes': signature_genes,
            'signature_scores': signature_scores,
            'n_programs': len(mp_programs),
            'patients_represented': list(set(p['patient'] for p in mp_programs))
        }
        print(f"    MP{mp_id} signature: {len(signature_genes)} genes, {len(meta_programs[mp_id]['patients_represented'])} samples")
        
    # Step 5: Calculate final Meta-Program scores for all cells
    print("\n📊 Step 5: Calculating final Meta-Program scores for all cells...")
    for mp_id, mp_data in meta_programs.items():
        available_signature_genes = [g for g in mp_data['signature_genes'] if g in malignant_cells.var_names]
        if len(available_signature_genes) >= 5:
            sc.tl.score_genes(malignant_cells, available_signature_genes, score_name=f'MP{mp_id}_score', use_raw=False)
            print(f"    MP{mp_id}: {len(available_signature_genes)} genes used for scoring")
    
    # Step 6: Assign cells to Meta-Programs
    print("\n🎯 Step 6: Assigning cells to Meta-Programs...")
    mp_score_cols = [f'MP{mp_id}_score' for mp_id in meta_programs.keys() if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(mp_score_cols) >= 2:
        score_df = malignant_cells.obs[mp_score_cols].copy()
        normalized_score_df = score_df.subtract(score_df.mean(axis=0), axis=1)
        
        mp_assignments = []
        for _, row in normalized_score_df.iterrows():
            sorted_scores = row.sort_values(ascending=False)
            max_score = sorted_scores.iloc[0]
            second_highest = sorted_scores.iloc[1]
            
            if second_highest <= 0.85 * max_score:
                assignment = sorted_scores.index[0].split('_')[0]
                mp_assignments.append(assignment)
            else:
                mp_assignments.append('Unresolved')
        
        malignant_cells.obs['MP_assignment'] = pd.Categorical(mp_assignments)
        
        print("    MP assignment distribution:")
        print(malignant_cells.obs['MP_assignment'].value_counts(normalize=True).mul(100).round(1).astype(str) + '%')

    # --- MODIFIED: Restored original save and return structure ---
    print("\n💾 Saving robust NMF results...")
    
    robust_nmf_results = {
        'all_programs': all_programs,
        'robust_programs': robust_programs,
        'meta_programs': meta_programs,
        'clustering_info': {
            'optimal_n_clusters': optimal_n_clusters,
            'similarity_matrix': similarity_matrix,
            'linkage_matrix': linkage_matrix,
            'cluster_assignments': clusters
        },
        'parameters': {
            'K_range': list(K_RANGE),
            'top_genes_per_program': TOP_GENES_PER_PROGRAM,
            'min_cells_per_patient': MIN_CELLS_PER_PATIENT,
            'n_meta_programs': N_META_PROGRAMS,
            'same_patient_overlap_threshold': SAME_PATIENT_OVERLAP_THRESHOLD,
            'cross_patient_overlap_threshold': CROSS_PATIENT_OVERLAP_THRESHOLD,
            'max_intra_patient_overlap': MAX_INTRA_PATIENT_OVERLAP
        }
    }
    
    # Save the results dictionary to a pickle file
    results_path = os.path.join(save_path, "robust_metabolism_nmf_results_aligned.pkl")
    with open(results_path, 'wb') as f:
        pickle.dump(robust_nmf_results, f)
    
    print(f"\n✅ Robust NMF analysis completed!")
    print(f"📊 Summary:")
    print(f"    - {len(all_programs)} total programs generated")
    print(f"    - {len(robust_programs)} robust programs identified")
    print(f"    - {optimal_n_clusters} Meta-Programs discovered")
    print(f"    - {len(valid_patients)} patients analyzed")
    
    return robust_nmf_results

# Manual Meta-Program Merger Function
def merge_meta_programs(malignant_cells, robust_results, 
                       source_mp=7, target_mp=1, 
                       save_path=None, update_scoring=True):
    """
    Manually merge one Meta-Program into another
    
    Parameters:
    -----------
    malignant_cells : AnnData
        Annotated data object with cells
    robust_results : dict
        Results from robust_nmf_analysis
    source_mp : int
        MP to be merged (will be eliminated), default 7
    target_mp : int
        MP to merge into (will be expanded), default 1
    save_path : str
        Path to save updated results
    update_scoring : bool
        Whether to recalculate MP scores after merger
    
    Returns:
    --------
    dict : Updated robust_results
    AnnData : Updated malignant_cells
    """
    
    print(f"🔄 Merging MP{source_mp} into MP{target_mp}...")
    print("="*50)
    
    # Create deep copy to avoid modifying original data
    updated_results = copy.deepcopy(robust_results)
    
    # =========================================================================
    # STEP 1: Update Meta-Programs Dictionary
    # =========================================================================
    
    print(f"📊 Step 1: Updating meta-programs dictionary...")
    
    meta_programs = updated_results['meta_programs']
    
    if source_mp not in meta_programs or target_mp not in meta_programs:
        print(f"❌ Error: MP{source_mp} or MP{target_mp} not found in meta_programs!")
        available_mps = list(meta_programs.keys())
        print(f"Available MPs: {available_mps}")
        return None, None
    
    # Get source and target MP data
    source_mp_data = meta_programs[source_mp]
    target_mp_data = meta_programs[target_mp]
    
    print(f"   MP{source_mp}: {len(source_mp_data['signature_genes'])} genes, "
          f"{source_mp_data['n_programs']} programs")
    print(f"   MP{target_mp}: {len(target_mp_data['signature_genes'])} genes, "
          f"{target_mp_data['n_programs']} programs")
    
    # Merge constituent programs
    merged_programs = target_mp_data['constituent_programs'] + source_mp_data['constituent_programs']
    
    # Merge signature genes with proper scoring
    merged_gene_data = merge_signature_genes(
        target_mp_data['signature_genes'], target_mp_data['signature_scores'],
        source_mp_data['signature_genes'], source_mp_data['signature_scores']
    )
    
    # Update target MP with merged data
    meta_programs[target_mp] = {
        'mp_id': target_mp,
        'constituent_programs': merged_programs,
        'signature_genes': merged_gene_data['genes'],
        'signature_scores': merged_gene_data['scores'],
        'n_programs': len(merged_programs),
        'patients_represented': list(set(
            target_mp_data['patients_represented'] + 
            source_mp_data['patients_represented']
        )),
        'merged_from': [source_mp],  # Track merger history
        'merge_timestamp': pd.Timestamp.now().isoformat()
    }
    
    # Remove source MP
    del meta_programs[source_mp]
    
    print(f"   ✅ Merged MP{target_mp}: {len(merged_gene_data['genes'])} genes, "
          f"{len(merged_programs)} programs, "
          f"{len(meta_programs[target_mp]['patients_represented'])} samples")
    
    # =========================================================================
    # STEP 2: Update Clustering Information
    # =========================================================================
    
    print(f"🌳 Step 2: Updating clustering information...")
    
    clustering_info = updated_results['clustering_info']
    
    # Update optimal number of clusters
    old_n_clusters = clustering_info['optimal_n_clusters']
    new_n_clusters = old_n_clusters - 1
    clustering_info['optimal_n_clusters'] = new_n_clusters
    
    # Update cluster assignments
    old_assignments = clustering_info['cluster_assignments'].copy()
    new_assignments = []
    
    for assignment in old_assignments:
        if assignment == source_mp:
            new_assignments.append(target_mp)  # Reassign source to target
        elif assignment > source_mp:
            new_assignments.append(assignment - 1)  # Shift down by 1
        else:
            new_assignments.append(assignment)  # Keep unchanged
    
    clustering_info['cluster_assignments'] = np.array(new_assignments)
    
    print(f"   ✅ Updated clusters: {old_n_clusters} → {new_n_clusters}")
    
    # =========================================================================
    # STEP 3: Renumber Meta-Programs (Optional but Recommended)
    # =========================================================================
    
    print(f"🔢 Step 3: Renumbering Meta-Programs...")
    
    # Create mapping for renumbering
    old_mp_ids = sorted([mp_id for mp_id in meta_programs.keys() if mp_id != source_mp])
    new_mp_mapping = {}
    
    new_id = 1
    for old_id in old_mp_ids:
        new_mp_mapping[old_id] = new_id
        new_id += 1
    
    # Renumber meta_programs dictionary
    new_meta_programs = {}
    for old_id, new_id in new_mp_mapping.items():
        mp_data = meta_programs[old_id].copy()
        mp_data['mp_id'] = new_id
        new_meta_programs[new_id] = mp_data
    
    updated_results['meta_programs'] = new_meta_programs
    
    print(f"   ✅ Renumbered MPs: {dict(new_mp_mapping)}")
    
    # =========================================================================
    # STEP 4: Update Cell Annotations
    # =========================================================================
    
    print(f"🧬 Step 4: Updating cell annotations...")
    
    # Update MP scores in malignant_cells
    if update_scoring:
        print("   Recalculating MP scores...")
        
        # Remove old MP scores
        old_score_cols = [col for col in malignant_cells.obs.columns 
                         if col.startswith('MP') and '_score' in col]
        for col in old_score_cols:
            del malignant_cells.obs[col]
        
        # Calculate new MP scores
        for new_id, mp_data in new_meta_programs.items():
            signature_genes = mp_data['signature_genes']
            available_genes = [g for g in signature_genes if g in malignant_cells.var_names]
            
            if len(available_genes) >= 5:
                sc.tl.score_genes(malignant_cells, available_genes, 
                                 score_name=f'MP{new_id}_score', use_raw=False)
                print(f"     MP{new_id}: {len(available_genes)} genes used for scoring")
    
    # Update MP assignments
    if 'MP_assignment' in malignant_cells.obs.columns:
        print("   Updating MP assignments...")
        
        old_assignments = malignant_cells.obs['MP_assignment'].copy()
        new_assignments = []
        
        for assignment in old_assignments:
            if assignment == f'MP{source_mp}':
                new_assignments.append(f'MP{new_mp_mapping[target_mp]}')
            elif assignment.startswith('MP'):
                try:
                    old_mp_num = int(assignment.replace('MP', ''))
                    if old_mp_num in new_mp_mapping:
                        new_assignments.append(f'MP{new_mp_mapping[old_mp_num]}')
                    else:
                        new_assignments.append('Unresolved')
                except:
                    new_assignments.append(assignment)  # Keep as is if can't parse
            else:
                new_assignments.append(assignment)  # Keep non-MP assignments
        
        malignant_cells.obs['MP_assignment'] = new_assignments
        malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
        
        # Print updated assignment statistics
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        print("   Updated MP assignment distribution:")
        for mp, count in assignment_counts.items():
            percentage = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} cells ({percentage:.1f}%)")
    
    # =========================================================================
    # STEP 5: Update Normalized Scores (if they exist)
    # =========================================================================
    
    print(f"📏 Step 5: Updating normalized scores...")
    
    # Remove old normalized scores
    normalized_cols = [col for col in malignant_cells.obs.columns 
                      if '_score_normalized' in col]
    for col in normalized_cols:
        del malignant_cells.obs[col]
    
    # Recalculate normalized scores
    new_score_cols = [f'MP{mp_id}_score' for mp_id in new_meta_programs.keys() 
                     if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(new_score_cols) >= 2:
        for col in new_score_cols:
            malignant_cells.obs[f'{col}_normalized'] = (
                malignant_cells.obs[col] - malignant_cells.obs[col].mean()
            )
        print(f"   ✅ Recalculated {len(new_score_cols)} normalized scores")
    
    # =========================================================================
    # STEP 6: Add Merger Metadata
    # =========================================================================
    
    print(f"📝 Step 6: Adding merger metadata...")
    
    # Add merger information to results
    if 'merger_history' not in updated_results:
        updated_results['merger_history'] = []
    
    merger_record = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'source_mp': source_mp,
        'target_mp': target_mp,
        'final_target_mp': new_mp_mapping[target_mp],
        'old_n_clusters': old_n_clusters,
        'new_n_clusters': new_n_clusters,
        'mp_mapping': new_mp_mapping,
        'merged_gene_count': len(merged_gene_data['genes']),
        'merged_program_count': len(merged_programs)
    }
    
    updated_results['merger_history'].append(merger_record)
    
    print(f"   ✅ Merger recorded in metadata")
    
    # =========================================================================
    # STEP 7: Save Updated Results
    # =========================================================================
    
    if save_path:
        print(f"💾 Step 7: Saving updated results...")
        
        # Save updated robust results
        with open(os.path.join(save_path, "robust_nmf_results_merged.pkl"), 'wb') as f:
            pickle.dump(updated_results, f)
        
        # Save updated malignant cells
        malignant_cells.write(os.path.join(save_path, "malignant_cells_merged_mps.h5ad"))
        
        # Save merger summary
        merger_summary = pd.DataFrame([merger_record])
        merger_summary.to_csv(os.path.join(save_path, "mp_merger_summary.csv"), index=False)
        
        print(f"   ✅ Results saved to {save_path}")
    
    # =========================================================================
    # STEP 8: Print Summary
    # =========================================================================
    
    print(f"\n🎉 Meta-Program Merger Completed!")
    print("="*50)
    print(f"✅ Merged MP{source_mp} into MP{target_mp}")
    print(f"📊 Meta-Programs: {old_n_clusters} → {new_n_clusters}")
    print(f"🧬 Final MP{new_mp_mapping[target_mp]} signature: {len(merged_gene_data['genes'])} genes")
    print(f"🔗 Merged programs: {len(merged_programs)} total")
    print(f"📍 Samples represented: {len(new_meta_programs[new_mp_mapping[target_mp]]['patients_represented'])}")
    
    return updated_results, malignant_cells

def merge_signature_genes(target_genes, target_scores, source_genes, source_scores, 
                         max_genes=50, score_weight=0.7):
    """
    Merge signature genes from two Meta-Programs with proper scoring
    
    Parameters:
    -----------
    target_genes, source_genes : list
        Gene lists from target and source MPs
    target_scores, source_scores : list  
        Corresponding gene scores
    max_genes : int
        Maximum number of genes in merged signature
    score_weight : float
        Weight for averaging scores (0.5 = equal weight)
    """
    
    # Create gene score dictionaries
    target_gene_scores = dict(zip(target_genes, target_scores))
    source_gene_scores = dict(zip(source_genes, source_scores))
    
    # Merge genes with weighted averaging for overlapping genes
    merged_gene_scores = {}
    
    # Add target genes
    for gene, score in target_gene_scores.items():
        merged_gene_scores[gene] = score
    
    # Add source genes (average scores for overlapping genes)
    for gene, score in source_gene_scores.items():
        if gene in merged_gene_scores:
            # Weighted average for overlapping genes
            existing_score = merged_gene_scores[gene]
            merged_gene_scores[gene] = (score_weight * existing_score + 
                                      (1 - score_weight) * score)
        else:
            # Add new gene from source
            merged_gene_scores[gene] = score
    
    # Sort by score and select top genes
    sorted_genes = sorted(merged_gene_scores.items(), key=lambda x: x[1], reverse=True)
    
    final_genes = [gene for gene, score in sorted_genes[:max_genes]]
    final_scores = [score for gene, score in sorted_genes[:max_genes]]
    
    return {
        'genes': final_genes,
        'scores': final_scores,
        'n_target_genes': len(target_genes),
        'n_source_genes': len(source_genes),
        'n_overlapping': len(set(target_genes) & set(source_genes)),
        'n_final_genes': len(final_genes)
    }

def validate_merger(updated_results, malignant_cells, save_path=None):
    """
    Validate the merger results and create summary statistics
    """
    
    print("🔍 Validating merger results...")
    
    validation_stats = {}
    
    # Check meta-programs
    meta_programs = updated_results['meta_programs']
    validation_stats['n_meta_programs'] = len(meta_programs)
    validation_stats['mp_ids'] = list(meta_programs.keys())
    
    # Check gene signature sizes
    signature_sizes = {mp_id: len(mp_data['signature_genes']) 
                      for mp_id, mp_data in meta_programs.items()}
    validation_stats['signature_sizes'] = signature_sizes
    
    # Check cell assignments
    if 'MP_assignment' in malignant_cells.obs.columns:
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        validation_stats['assignment_counts'] = assignment_counts.to_dict()
        
        # Check for unresolved cells
        unresolved_pct = (assignment_counts.get('Unresolved', 0) / 
                         len(malignant_cells) * 100)
        validation_stats['unresolved_percentage'] = unresolved_pct
    
    # Check MP scores
    mp_score_cols = [col for col in malignant_cells.obs.columns 
                    if col.startswith('MP') and '_score' in col and 'normalized' not in col]
    validation_stats['available_mp_scores'] = mp_score_cols
    
    # Print validation summary
    print("📊 Validation Summary:")
    print(f"   Meta-Programs: {validation_stats['n_meta_programs']}")
    print(f"   MP IDs: {validation_stats['mp_ids']}")
    print(f"   Signature sizes: {signature_sizes}")
    print(f"   Available MP scores: {len(mp_score_cols)}")
    
    if 'assignment_counts' in validation_stats:
        print(f"   Cell assignments:")
        for mp, count in validation_stats['assignment_counts'].items():
            pct = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} ({pct:.1f}%)")
        print(f"   Unresolved: {unresolved_pct:.1f}%")
    
    if save_path:
        validation_df = pd.DataFrame([validation_stats])
        validation_df.to_csv(os.path.join(save_path, "merger_validation.csv"), index=False)
    
    return validation_stats

def split_meta_program(malignant_cells, robust_results, 
                      mp_to_split, n_splits=2, 
                      split_method='expression_clustering',
                      save_path=None, update_scoring=True,
                      min_cells_per_split=50, min_genes_per_split=10):
    """
    Split a Meta-Program into multiple distinct programs
    
    Parameters:
    -----------
    malignant_cells : AnnData
        Annotated data object with cells
    robust_results : dict
        Results from robust_nmf_analysis
    mp_to_split : int
        MP to be split
    n_splits : int
        Number of programs to split into
    split_method : str
        Method for splitting ('expression_clustering', 'program_clustering', 'gene_clustering')
    save_path : str
        Path to save updated results
    update_scoring : bool
        Whether to recalculate MP scores after split
    min_cells_per_split : int
        Minimum cells required per split
    min_genes_per_split : int
        Minimum genes required per split
    
    Returns:
    --------
    dict : Updated robust_results
    AnnData : Updated malignant_cells
    """
    
    print(f"✂️ Splitting MP{mp_to_split} into {n_splits} programs...")
    print("="*60)
    
    # Create deep copy to avoid modifying original data
    updated_results = copy.deepcopy(robust_results)
    
    # =========================================================================
    # STEP 1: Validate Input and Extract MP Data
    # =========================================================================
    
    print(f"📊 Step 1: Validating input and extracting MP data...")
    
    meta_programs = updated_results['meta_programs']
    
    if mp_to_split not in meta_programs:
        print(f"❌ Error: MP{mp_to_split} not found in meta_programs!")
        available_mps = list(meta_programs.keys())
        print(f"Available MPs: {available_mps}")
        return None, None
    
    mp_data = meta_programs[mp_to_split]
    
    print(f"   MP{mp_to_split}: {len(mp_data['signature_genes'])} genes, "
          f"{mp_data['n_programs']} programs, "
          f"{len(mp_data['patients_represented'])} samples")
    
    # Check if MP has enough data to split
    if mp_data['n_programs'] < n_splits:
        print(f"❌ Error: MP{mp_to_split} has only {mp_data['n_programs']} programs, "
              f"cannot split into {n_splits}")
        return None, None
    
    # Get cells assigned to this MP
    mp_cells_mask = malignant_cells.obs['MP_assignment'] == f'MP{mp_to_split}'
    n_mp_cells = mp_cells_mask.sum()
    
    print(f"   Cells assigned to MP{mp_to_split}: {n_mp_cells}")
    
    if n_mp_cells < min_cells_per_split * n_splits:
        print(f"❌ Error: Not enough cells for splitting. Need {min_cells_per_split * n_splits}, "
              f"have {n_mp_cells}")
        return None, None
    
    # =========================================================================
    # STEP 2: Perform Splitting Based on Selected Method
    # =========================================================================
    
    print(f"✂️ Step 2: Performing split using '{split_method}' method...")
    
    if split_method == 'program_clustering':
        split_results = split_by_program_clustering(
            mp_data, n_splits, min_genes_per_split
        )
    elif split_method == 'gene_clustering':
        split_results = split_by_gene_clustering(
            malignant_cells, mp_data, mp_cells_mask, n_splits, 
            min_genes_per_split, save_path
        )
    else:
        print(f"❌ Error: Unknown split method '{split_method}'")
        return None, None
    
    if split_results is None:
        print(f"❌ Splitting failed with method '{split_method}'")
        return None, None
    
    # =========================================================================
    # STEP 3: Create New Meta-Programs from Split Results
    # =========================================================================
    
    print(f"🧬 Step 3: Creating new Meta-Programs from split results...")
    
    # Remove original MP
    del meta_programs[mp_to_split]
    
    # Find the next available MP ID
    existing_mp_ids = list(meta_programs.keys())
    next_mp_id = max(existing_mp_ids) + 1 if existing_mp_ids else 1
    
    new_mp_data = {}
    new_mp_mapping = {}  # Track old MP -> new MPs
    
    for i, split_data in enumerate(split_results):
        new_mp_id = next_mp_id + i
        
        # Create new MP data structure
        new_mp_data[new_mp_id] = {
            'mp_id': new_mp_id,
            'constituent_programs': split_data['programs'],
            'signature_genes': split_data['signature_genes'],
            'signature_scores': split_data['signature_scores'],
            'n_programs': len(split_data['programs']),
            'patients_represented': split_data['patients_represented'],
            'split_from': mp_to_split,
            'split_method': split_method,
            'split_index': i,
            'split_timestamp': pd.Timestamp.now().isoformat()
        }
        
        print(f"   MP{new_mp_id}: {len(split_data['signature_genes'])} genes, "
              f"{len(split_data['programs'])} programs, "
              f"{len(split_data['patients_represented'])} samples")
    
    # Add new MPs to meta_programs
    meta_programs.update(new_mp_data)
    new_mp_mapping[mp_to_split] = list(new_mp_data.keys())
    
    # =========================================================================
    # STEP 4: Update Clustering Information
    # =========================================================================
    
    print(f"🌳 Step 4: Updating clustering information...")
    
    clustering_info = updated_results['clustering_info']
    
    # Update number of clusters
    old_n_clusters = clustering_info['optimal_n_clusters']
    new_n_clusters = old_n_clusters + n_splits - 1
    clustering_info['optimal_n_clusters'] = new_n_clusters
    
    # Update cluster assignments (this is complex for splits)
    old_assignments = clustering_info['cluster_assignments'].copy()
    
    # For simplicity, we'll mark split assignments as requiring recalculation
    clustering_info['cluster_assignments_valid'] = False
    clustering_info['requires_reclustering'] = True
    
    print(f"   ✅ Updated clusters: {old_n_clusters} → {new_n_clusters}")
    print(f"   ⚠️ Cluster assignments marked for recalculation")
    
    # =========================================================================
    # STEP 5: Update Cell Annotations
    # =========================================================================
    
    print(f"🧬 Step 5: Updating cell annotations...")
    
    # Update MP scores
    if update_scoring:
        print("   Recalculating MP scores...")
        
        # Remove old MP scores for split MP
        old_score_cols = [col for col in malignant_cells.obs.columns 
                         if col.startswith(f'MP{mp_to_split}_score')]
        for col in old_score_cols:
            del malignant_cells.obs[col]
        
        # Calculate new MP scores for split MPs
        for new_mp_id, mp_data in new_mp_data.items():
            signature_genes = mp_data['signature_genes']
            available_genes = [g for g in signature_genes if g in malignant_cells.var_names]
            
            if len(available_genes) >= min_genes_per_split:
                sc.tl.score_genes(malignant_cells, available_genes, 
                                 score_name=f'MP{new_mp_id}_score', use_raw=False)
                print(f"     MP{new_mp_id}: {len(available_genes)} genes used for scoring")
    
    # Update MP assignments using new scores
    if 'cell_assignments' in split_results[0]:
        print("   Updating MP assignments based on split results...")
        
        # Update assignments for cells that were in the split MP
        mp_cells_indices = np.where(mp_cells_mask)[0]
        
        for i, split_data in enumerate(split_results):
            new_mp_id = list(new_mp_data.keys())[i]
            
            # Get cell indices for this split
            if 'cell_assignments' in split_data:
                split_cell_indices = split_data['cell_assignments']
                
                # Update assignments
                for cell_idx in split_cell_indices:
                    if cell_idx < len(malignant_cells):
                        malignant_cells.obs.iloc[cell_idx, 
                            malignant_cells.obs.columns.get_loc('MP_assignment')] = f'MP{new_mp_id}'
    
    else:
        # Reassign cells based on highest MP scores
        print("   Reassigning cells based on highest scores...")
        reassign_cells_after_split(malignant_cells, mp_to_split, list(new_mp_data.keys()))
    
    # Update assignment category
    malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].astype('category')
    
    # Print updated assignment statistics
    assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
    print("   Updated MP assignment distribution:")
    for mp, count in assignment_counts.items():
        if mp.startswith('MP') and any(str(mp_id) in mp for mp_id in new_mp_data.keys()):
            percentage = count / len(malignant_cells) * 100
            print(f"     {mp}: {count} cells ({percentage:.1f}%)")
    
    # =========================================================================
    # STEP 6: Add Split Metadata
    # =========================================================================
    
    print(f"📝 Step 6: Adding split metadata...")
    
    # Add split information to results
    if 'split_history' not in updated_results:
        updated_results['split_history'] = []
    
    split_record = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'original_mp': mp_to_split,
        'split_method': split_method,
        'n_splits': n_splits,
        'new_mp_ids': list(new_mp_data.keys()),
        'old_n_clusters': old_n_clusters,
        'new_n_clusters': new_n_clusters,
        'split_gene_counts': [len(mp['signature_genes']) for mp in new_mp_data.values()],
        'split_program_counts': [mp['n_programs'] for mp in new_mp_data.values()]
    }
    
    updated_results['split_history'].append(split_record)
    
    print(f"   ✅ Split recorded in metadata")
    
    # =========================================================================
    # STEP 7: Save Updated Results
    # =========================================================================
    
    if save_path:
        print(f"💾 Step 7: Saving updated results...")
        
        # Save updated robust results
        with open(os.path.join(save_path, "robust_nmf_results_split.pkl"), 'wb') as f:
            pickle.dump(updated_results, f)
        
        # Save updated malignant cells
        malignant_cells.write(os.path.join(save_path, "malignant_cells_split_mps.h5ad"))
        
        # Save split summary
        split_summary = pd.DataFrame([split_record])
        split_summary.to_csv(os.path.join(save_path, "mp_split_summary.csv"), index=False)
        
        print(f"   ✅ Results saved to {save_path}")
    
    # =========================================================================
    # STEP 8: Print Summary
    # =========================================================================
    
    print(f"\n🎉 Meta-Program Split Completed!")
    print("="*60)
    print(f"✂️ Split MP{mp_to_split} into {n_splits} new MPs")
    print(f"📊 Meta-Programs: {old_n_clusters} → {new_n_clusters}")
    print(f"🆕 New MP IDs: {list(new_mp_data.keys())}")
    
    for new_mp_id, mp_data in new_mp_data.items():
        print(f"   MP{new_mp_id}: {len(mp_data['signature_genes'])} genes, "
              f"{mp_data['n_programs']} programs")
    
    return updated_results, malignant_cells

def split_by_program_clustering(mp_data, n_splits, min_genes_per_split):
    """
    Split MP based on clustering of constituent programs
    """
    
    print(f"   🔧 Using program-based clustering...")
    
    constituent_programs = mp_data['constituent_programs']
    
    if len(constituent_programs) < n_splits:
        print(f"   ❌ Not enough constituent programs for clustering")
        return None
    
    # Create similarity matrix between programs
    n_programs = len(constituent_programs)
    similarity_matrix = np.zeros((n_programs, n_programs))
    
    for i in range(n_programs):
        for j in range(n_programs):
            genes_i = set(constituent_programs[i]['genes'])
            genes_j = set(constituent_programs[j]['genes'])
            jaccard_sim = len(genes_i & genes_j) / len(genes_i | genes_j) if genes_i | genes_j else 0
            similarity_matrix[i, j] = jaccard_sim
    
    # Hierarchical clustering
    distance_matrix = 1 - similarity_matrix
    condensed_distances = pdist(distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='ward')
    cluster_labels = fcluster(linkage_matrix, n_splits, criterion='maxclust')
    
    # Create split results
    split_results = []
    
    for cluster_id in range(1, n_splits + 1):
        cluster_programs = [constituent_programs[i] for i, label in enumerate(cluster_labels) 
                           if label == cluster_id]
        
        if len(cluster_programs) == 0:
            continue
        
        # Merge genes from cluster programs
        all_genes = []
        all_scores = []
        
        for program in cluster_programs:
            all_genes.extend(program['genes'])
            all_scores.extend(program.get('gene_scores', [1.0] * len(program['genes'])))
        
        # Count gene frequencies and calculate average scores
        gene_counts = Counter(all_genes)
        gene_score_sums = {}
        
        for gene, score in zip(all_genes, all_scores):
            if gene not in gene_score_sums:
                gene_score_sums[gene] = []
            gene_score_sums[gene].append(score)
        
        # Select genes appearing in multiple programs
        cluster_genes = []
        cluster_scores = []
        
        for gene, count in gene_counts.items():
            if count >= max(1, len(cluster_programs) * 0.3):  # 30% threshold
                avg_score = np.mean(gene_score_sums[gene])
                cluster_genes.append(gene)
                cluster_scores.append(avg_score)
        
        # Sort by score and take top genes
        if len(cluster_genes) >= min_genes_per_split:
            gene_score_pairs = list(zip(cluster_genes, cluster_scores))
            gene_score_pairs.sort(key=lambda x: x[1], reverse=True)
            
            final_genes = [gene for gene, score in gene_score_pairs[:40]]
            final_scores = [score for gene, score in gene_score_pairs[:40]]
            
            split_results.append({
                'signature_genes': final_genes,
                'signature_scores': final_scores,
                'programs': cluster_programs,
                'patients_represented': list(set(p['patient'] for p in cluster_programs)),
                'n_programs': len(cluster_programs)
            })
    
    return split_results

def split_by_gene_clustering(malignant_cells, mp_data, mp_cells_mask, 
                           n_splits, min_genes_per_split, save_path=None):
    """
    Split MP based on clustering of signature genes
    """
    
    print(f"   🧬 Using gene-based clustering...")
    
    signature_genes = mp_data['signature_genes']
    
    if len(signature_genes) < min_genes_per_split * n_splits:
        print(f"   ❌ Not enough genes for clustering")
        return None
    
    # Get expression of signature genes in MP cells
    mp_cells = malignant_cells[mp_cells_mask]
    available_genes = [g for g in signature_genes if g in mp_cells.var_names]
    
    if len(available_genes) < min_genes_per_split * n_splits:
        print(f"   ❌ Not enough available genes for clustering")
        return None
    
    # Calculate gene-gene correlation matrix
    X_genes = mp_cells[:, available_genes].X.T  # Genes x Cells
    if hasattr(X_genes, 'toarray'):
        X_genes = X_genes.toarray()
    
    # Compute pairwise correlations
    gene_corr_matrix = np.corrcoef(X_genes)
    gene_corr_matrix = np.nan_to_num(gene_corr_matrix)
    
    # Convert correlation to distance
    gene_distance_matrix = 1 - np.abs(gene_corr_matrix)
    
    # Hierarchical clustering of genes
    condensed_distances = pdist(gene_distance_matrix)
    linkage_matrix = linkage(condensed_distances, method='ward')
    gene_cluster_labels = fcluster(linkage_matrix, n_splits, criterion='maxclust')
    
    # Create split results based on gene clusters
    split_results = []
    
    for cluster_id in range(1, n_splits + 1):
        cluster_gene_indices = np.where(gene_cluster_labels == cluster_id)[0]
        cluster_genes = [available_genes[i] for i in cluster_gene_indices]
        
        if len(cluster_genes) < min_genes_per_split:
            continue
        
        # Calculate scores for genes in this cluster
        cluster_scores = [mp_data['signature_scores'][signature_genes.index(gene)] 
                         for gene in cluster_genes if gene in signature_genes]
        
        # Pad scores if needed
        while len(cluster_scores) < len(cluster_genes):
            cluster_scores.append(1.0)
        
        # Distribute constituent programs
        n_programs_for_cluster = len(mp_data['constituent_programs']) // n_splits
        start_idx = (cluster_id - 1) * n_programs_for_cluster
        end_idx = start_idx + n_programs_for_cluster if cluster_id < n_splits else len(mp_data['constituent_programs'])
        cluster_programs = mp_data['constituent_programs'][start_idx:end_idx]
        
        split_results.append({
            'signature_genes': cluster_genes,
            'signature_scores': cluster_scores,
            'programs': cluster_programs,
            'patients_represented': list(set(p['patient'] for p in cluster_programs)),
            'gene_cluster_id': cluster_id
        })
    
    return split_results

def reassign_cells_after_split(malignant_cells, original_mp, new_mp_ids):
    """
    Reassign cells originally assigned to split MP based on new MP scores
    Fixed version that handles Categorical columns properly
    """
    
    # Get cells that were assigned to original MP
    original_assignment_mask = malignant_cells.obs['MP_assignment'] == f'MP{original_mp}'
    
    if original_assignment_mask.sum() == 0:
        print(f"   ⚠️ No cells found assigned to MP{original_mp}")
        return
    
    print(f"   📊 Found {original_assignment_mask.sum()} cells assigned to MP{original_mp}")
    
    # Get scores for new MPs
    new_score_cols = [f'MP{mp_id}_score' for mp_id in new_mp_ids 
                     if f'MP{mp_id}_score' in malignant_cells.obs.columns]
    
    if len(new_score_cols) == 0:
        print("   ⚠️ No new MP scores available for reassignment")
        return
    
    print(f"   🔄 Using scores from: {new_score_cols}")
    
    # SOLUTION 1: Convert to object type temporarily, then back to categorical
    # This is the most robust approach
    
    # Store original categorical info
    original_categories = malignant_cells.obs['MP_assignment'].cat.categories.tolist()
    
    # Add new MP categories if they don't exist
    new_categories = [f'MP{mp_id}' for mp_id in new_mp_ids]
    all_categories = list(set(original_categories + new_categories))
    
    # Update categories first
    malignant_cells.obs['MP_assignment'] = malignant_cells.obs['MP_assignment'].cat.add_categories(
        [cat for cat in new_categories if cat not in original_categories]
    )
    
    # Now reassign cells based on highest score
    reassignment_counts = {f'MP{mp_id}': 0 for mp_id in new_mp_ids}
    
    for cell_idx in np.where(original_assignment_mask)[0]:
        # Get scores for this cell
        scores = []
        for col in new_score_cols:
            score = malignant_cells.obs.iloc[cell_idx][col]
            scores.append(score)
        
        # Find MP with highest score
        best_mp_idx = np.argmax(scores)
        best_mp_id = new_mp_ids[best_mp_idx]
        best_score = scores[best_mp_idx]
        
        # Assign cell to best MP
        malignant_cells.obs.iloc[cell_idx, 
            malignant_cells.obs.columns.get_loc('MP_assignment')] = f'MP{best_mp_id}'
        
        reassignment_counts[f'MP{best_mp_id}'] += 1
    
    # Print reassignment summary
    print(f"   ✅ Reassignment completed:")
    for mp_name, count in reassignment_counts.items():
        if count > 0:
            percentage = count / original_assignment_mask.sum() * 100
            print(f"     {mp_name}: {count} cells ({percentage:.1f}%)")

def create_split_visualization(X_data, cluster_labels, split_results, title, save_path):
    """
    Create visualization of the split results
    """
    
    # PCA for visualization
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_data)
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: PCA with cluster colors
    scatter = axes[0].scatter(X_pca[:, 0], X_pca[:, 1], c=cluster_labels, 
                             cmap='tab10', alpha=0.7)
    axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[0].set_title(f'{title} - Cell Clustering')
    plt.colorbar(scatter, ax=axes[0])
    
    # Plot 2: Split statistics
    split_names = [f'Split {i+1}' for i in range(len(split_results))]
    gene_counts = [len(result['signature_genes']) for result in split_results]
    
    bars = axes[1].bar(split_names, gene_counts, alpha=0.7)
    axes[1].set_ylabel('Number of Signature Genes')
    axes[1].set_title('Split Results - Gene Counts')
    
    # Add value labels on bars
    for bar, count in zip(bars, gene_counts):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    str(count), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f'{save_path}/{title}_split_visualization.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{save_path}/{title}_split_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()

def export_results_for_r_analysis(malignant_cells, robust_results, save_path):
    """
    Export all necessary results for R statistical analysis and visualization
    """
 
    # Create R analysis directory
    r_analysis_dir = os.path.join(save_path, "NMF_results_for_R")
    os.makedirs(r_analysis_dir, exist_ok=True)
    
    # 1. Clinical metadata with MP assignments
    print("📊 Exporting clinical metadata...")
    clinical_data = malignant_cells.obs.copy()
    
    # Add MP scores if available
    mp_score_cols = [col for col in clinical_data.columns if 'MP' in col and 'score' in col]
    
    # Ensure all required columns are present
    required_cols = ['patient', 'sample_id', 'tissue', 'study']
    
    # Add MP assignment columns
    mp_cols = [col for col in clinical_data.columns if col.startswith('MP') or 'Meta' in col]
    
    export_cols = [col for col in required_cols + mp_cols + mp_score_cols if col in clinical_data.columns]
    
    clinical_export = clinical_data[export_cols].copy()
    clinical_export.to_csv(os.path.join(r_analysis_dir, "clinical_metadata.csv"), index=True)
    
    # 2. Patient-level summary statistics
    print("📈 Creating patient-level summary...")
    patient_summary = []
    
    for patient_id in clinical_data['patient'].unique():
        patient_cells = clinical_data[clinical_data['patient'] == patient_id]
        
        # Basic info
        patient_info = {
            'patient': patient_id,
            'n_cells': len(patient_cells),
            'study': patient_cells['study'].iloc[0] if 'study' in patient_cells.columns else 'Unknown'
        }
        
        # MP proportions if available
        if 'MP_assignment' in patient_cells.columns:
            mp_counts = patient_cells['MP_assignment'].value_counts()
            total_cells = len(patient_cells)
            
            for mp in mp_counts.index:
                patient_info[f'{mp}_count'] = mp_counts[mp]
                patient_info[f'{mp}_proportion'] = mp_counts[mp] / total_cells
        
        # MP scores if available
        for col in mp_score_cols:
            if col in patient_cells.columns:
                patient_info[f'{col}_mean'] = patient_cells[col].mean()
        
        patient_summary.append(patient_info)
    
    patient_df = pd.DataFrame(patient_summary)
    patient_df.to_csv(os.path.join(r_analysis_dir, "patient_level_summary.csv"), index=False)
    
    # 3. MP signature genes if available
    print("🧬 Exporting MP signature genes...")
    if robust_results and 'meta_programs' in robust_results:
        mp_signatures = robust_results['meta_programs']
        
        # Create signature matrix
        all_genes = set()
        for mp_name, genes in mp_signatures.items():
            all_genes.update(genes)
        
        signature_matrix = pd.DataFrame(0, index=list(all_genes), columns=list(mp_signatures.keys()))
        
        for mp_name, genes in mp_signatures.items():
            for gene in genes:
                signature_matrix.loc[gene, mp_name] = 1
        
        signature_matrix.to_csv(os.path.join(r_analysis_dir, "mp_signature_genes.csv"))
        
        # Export individual signatures as separate files
        for mp_name, genes in mp_signatures.items():
            pd.DataFrame({'gene': genes}).to_csv(
                os.path.join(r_analysis_dir, f"{mp_name}_signature_genes.csv"), index=False)
    
    # 4. Save expression matrix in efficient sparse format (like scanpy standard)
    print("🔬 Exporting expression data in sparse matrix format...")
    
    # Get key genes (MP signatures + highly variable genes)
    key_genes = set()
    
    # Add MP signature genes
    if robust_results and 'meta_programs' in robust_results:
        for program_ in robust_results['meta_programs'].keys():
            key_genes.update(robust_results['meta_programs'][program_]['signature_genes'])
    
    # Add highly variable genes if available
    if hasattr(malignant_cells, 'var') and 'highly_variable' in malignant_cells.var.columns:
        hvg_genes = malignant_cells.var[malignant_cells.var['highly_variable']].index.tolist()
        key_genes.update(hvg_genes[:1000])  # Top 1000 HVGs
    
    # Filter genes that exist in the data
    available_genes = [g for g in key_genes if g in malignant_cells.var.index]
    
    print(f"   Exporting {len(available_genes)} genes in sparse format...")
    
    if available_genes:
        # Subset AnnData object to key genes
        adata_subset = malignant_cells[:, available_genes].copy()
        
        # Convert expression matrix to sparse CSR format and save as .mtx
        # Ensure matrix is in sparse format
        if scipy.sparse.issparse(adata_subset.X):
            matrix = scipy.sparse.csr_matrix(adata_subset.X)
        else:
            matrix = scipy.sparse.csr_matrix(np.array(adata_subset.X))
        
        # Save expression matrix in Matrix Market format (.mtx)
        scipy.io.mmwrite(os.path.join(r_analysis_dir, "expression_matrix.mtx"), matrix)
        
        # Save gene names (rows)
        adata_subset.var_names.to_series().to_csv(
            os.path.join(r_analysis_dir, "gene_names.csv"), 
            header=['gene_name'], index=False
        )
        
        # Save cell barcodes (columns) 
        adata_subset.obs_names.to_series().to_csv(
            os.path.join(r_analysis_dir, "cell_barcodes.csv"), 
            header=['barcode'], index=False
        )
        
        # Save additional gene metadata if available
        if adata_subset.var.shape[1] > 0:
            adata_subset.var.to_csv(os.path.join(r_analysis_dir, "gene_metadata.csv"))
        
        print(f"   ✅ Sparse matrix saved: {matrix.shape[0]} cells × {matrix.shape[1]} genes")
        print(f"   💾 Matrix sparsity: {(1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])):.1%}")
    
    # 5. NMF analysis parameters and results summary
    print("📋 Exporting analysis parameters...")
    if robust_results:
        analysis_summary = {
            'n_cells': malignant_cells.n_obs,
            'n_genes': malignant_cells.n_vars,
            'n_patients': clinical_data['patient'].nunique(),
            'n_meta_programs': len(robust_results.get('meta_programs', {})),
        }
        
        # Add parameter information
        if 'parameters' in robust_results:
            analysis_summary.update(robust_results['parameters'])
        
        pd.DataFrame([analysis_summary]).to_csv(
            os.path.join(r_analysis_dir, "analysis_summary.csv"), index=False)
        
    # 6. Create data dictionary
    print("📚 Creating data dictionary...")
    data_dict = {
        'clinical_metadata.csv': 'Cell-level clinical metadata with MP assignments and scores',
        'patient_level_summary.csv': 'Patient-level aggregated data with MP proportions and statistics',
        'mp_signature_genes.csv': 'Binary matrix of genes in each meta-program signature',
        'expression_matrix.mtx': 'Sparse expression matrix in Matrix Market format (cells × genes)',
        'gene_names.csv': 'Gene names corresponding to expression matrix rows',
        'cell_barcodes.csv': 'Cell barcodes corresponding to expression matrix columns',
        'gene_metadata.csv': 'Additional gene metadata (if available)',
        'malignant_cells_complete.h5ad': 'Complete AnnData object with all data',
        'analysis_summary.csv': 'Summary of analysis parameters and basic statistics',
        'msi_treated_cells.csv': 'Subset of MSI treated cells for response analysis',
        'msi_treated_patients.csv': 'Patient-level data for MSI treated patients only'
    }
    
    pd.DataFrame(list(data_dict.items()), columns=['File', 'Description']).to_csv(
        os.path.join(r_analysis_dir, "data_dictionary.csv"), index=False)
    
    # NEW EXPORT 1: MP fraction data for stacked barplot
    print("📊 Exporting MP fraction data for stacked barplot...")
    if 'MP_assignment' in malignant_cells.obs.columns:
        sample_mp_fractions = calculate_mp_fractions_per_sample(malignant_cells, 'MP_assignment')
        sample_mp_fractions.to_csv(os.path.join(r_analysis_dir, "sample_mp_fractions.csv"), index=False)
        print(f"   ✅ MP fractions: {sample_mp_fractions.shape}")
    
    # NEW EXPORT 2: MP marker genes with rankings for heatmap
    print("🧬 Exporting MP marker genes for heatmap...")
    if robust_results and 'meta_programs' in robust_results:
        mp_marker_data = []
        for mp_name, program_ in robust_results['meta_programs'].items():
            for i, gene in enumerate(program_["signature_genes"][:50]):  # Top 50 genes per MP
                mp_marker_data.append({
                    'gene': gene,
                    'Meta_Program': mp_name,
                    'rank': i + 1,
                    'weight': 1.0 - (i * 0.02)  # Decreasing importance
                })
        
        if mp_marker_data:
            mp_marker_df = pd.DataFrame(mp_marker_data)
            mp_marker_df.to_csv(os.path.join(r_analysis_dir, "mp_marker_genes_ranked.csv"), index=False)
            print(f"   ✅ MP marker genes: {len(mp_marker_data)} gene-MP pairs")
    
    # NEW EXPORT 3: Expression subset for heatmap (top 20 genes per MP)
    print("🔥 Exporting expression subset for marker gene heatmap...")
    if robust_results and 'meta_programs' in robust_results:
        # Get top 20 genes per MP
        top_genes = set()
        for program_ in robust_results['meta_programs'].values():
            top_genes.update(program_["signature_genes"][:20])
        
        top_genes_list = list(top_genes)
        
        # Check which genes are available in expression matrix
        available_genes = []
        gene_indices = []
        
        for i, gene in enumerate(malignant_cells.var_names):
            if gene in top_genes_list:
                available_genes.append(gene)
                gene_indices.append(i)
        
        if available_genes:
            # Extract expression subset
            expr_subset = malignant_cells.X[:, gene_indices]
            
            # Convert to dense if sparse
            if hasattr(expr_subset, 'toarray'):
                expr_subset = expr_subset.toarray()
            
            # Create DataFrame (genes as rows, cells as columns)
            expr_subset_df = pd.DataFrame(
                expr_subset.T,  # Transpose so genes are rows
                index=available_genes,
                columns=malignant_cells.obs_names
            )
            
            expr_subset_df.to_csv(os.path.join(r_analysis_dir, "expression_matrix_subset.csv"))
            print(f"   ✅ Expression subset: {expr_subset_df.shape} (genes × cells)")
    
    # Update data dictionary
    data_dict.update({
        'sample_mp_fractions.csv': 'MP fractions per sample for stacked barplot',
        'mp_marker_genes_ranked.csv': 'Ranked marker genes per MP for heatmap',
        'expression_matrix_subset.csv': 'Expression data for top MP marker genes'
    })
    
    print(f"\n✅ Enhanced data export completed!")
    print(f"📁 Total files exported: {len(data_dict)}")
    
    return r_analysis_dir

# UMAP Visualization of MP Distribution
def create_mp_umap_visualization(malignant_cells, available_metabolism_genes, save_path, figsize=(15, 10)):
    """
    Create UMAP visualization showing MP distribution across cells.
    This corrected version improves subplot management to ensure figures are placed correctly
    by using plt.subplots() for better axis control and closing figures after use.

    Parameters:
    -----------
    malignant_cells : AnnData
        Malignant cells with MP assignments and scores.
    available_metabolism_genes : list
        List of metabolism genes to use for UMAP.
    save_path : str
        Directory to save figures.
    figsize : tuple
        Default figure size, can be overridden by specific plot settings.
    """
    print("🗺️ Creating UMAP visualization of Meta-Program distribution...")

    # Create a copy for UMAP analysis
    adata_umap = malignant_cells.copy()
    adata_umap = adata_umap[:, available_metabolism_genes]

    # Downsample if too many cells (for computational efficiency)
    if adata_umap.n_obs > 30000:
        print(f"  📉 Downsampling from {adata_umap.n_obs} to 30,000 cells for UMAP...")
        sc.pp.subsample(adata_umap, n_obs=30000, random_state=42)

    # Preprocessing for UMAP
    print("  🔄 Preprocessing for UMAP...")
    if adata_umap.X.max() > 50:  # Check if data looks like raw counts
        sc.pp.normalize_total(adata_umap, target_sum=1e4)
        sc.pp.log1p(adata_umap)
    sc.pp.scale(adata_umap, max_value=10)
    print("  📊 Computing PCA...")
    sc.tl.pca(adata_umap, svd_solver='arpack', n_comps=50)
    print("  🕸️ Computing neighborhood graph...")
    sc.pp.neighbors(adata_umap, n_neighbors=15, n_pcs=40)
    print("  🗺️ Computing UMAP embedding...")
    sc.tl.umap(adata_umap, random_state=42)

    # --- Plotting Section ---
    print("  🎨 Creating UMAP visualizations...")
    sc.settings.set_figure_params(dpi=300, facecolor='white')

    # --- Plot 1: Comprehensive UMAP overview (MP Assignment, Patient, Treatment) ---
    # CORRECTION: Use plt.subplots for better control over the figure and axes.
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    
    plot_configs = [
        {'key': 'MP_assignment', 'title': 'Meta-Program Assignment', 'palette': 'tab20', 'size': 30, 'alpha': 0.8, 'legend_loc': 'on data'},
        {'key': 'patient', 'title': 'Patient ID', 'palette': 'plasma', 'size': 20, 'alpha': 0.6, 'legend_loc': None},
        {'key': 'Treatment_Strategy', 'title': 'Treatment Strategy', 'palette': 'Set2', 'size': 30, 'alpha': 0.8, 'legend_loc': 'on data'}
    ]

    for i, config in enumerate(plot_configs):
        if config['key'] in adata_umap.obs.columns:
            sc.pl.umap(adata_umap, color=config['key'],
                       palette=config['palette'], size=config['size'], alpha=config['alpha'],
                       title=config['title'],
                       frameon=False, ax=axes[i], show=False, legend_loc=config['legend_loc'])
        else:
            axes[i].axis('off')
            axes[i].text(0.5, 0.5, f"'{config['key']}' not found", ha='center', va='center', transform=axes[i].transAxes)

    plt.tight_layout()
    plt.savefig(f"{save_path}/mp_umap_overview.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig) # CORRECTION: Close the figure to free memory and prevent state leakage.

    # --- Plot 2: UMAP colored by MP scores ---
    mp_score_cols = [col for col in adata_umap.obs.columns if 'MP' in col and 'score' in col and not 'normalized' in col]

    if mp_score_cols:
        all_scores = adata_umap.obs[mp_score_cols].values.flatten()
        vmin, vmax = np.nanmin(all_scores), np.nanmax(all_scores)

        # CORRECTION: Robustly calculate grid size to fit all plots.
        n_plots = len(mp_score_cols)
        ncols = 3 
        nrows = int(np.ceil(n_plots / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)
        axes = axes.flatten()

        for i, mp_col in enumerate(mp_score_cols):
            sc.pl.umap(adata_umap, color=mp_col,
                       size=5, alpha=0.8, cmap='viridis', # Using perceptually uniform colormap
                       title=f'{mp_col.replace("_", " ").title()}',
                       frameon=False, ax=axes[i], show=False,
                       vmin=vmin, vmax=vmax)

        # CORRECTION: Hide any unused subplots for a cleaner look.
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        plt.tight_layout()
        plt.savefig(f"{save_path}/mp_scores_umap.pdf", dpi=300, bbox_inches='tight')
        plt.show()
        plt.close(fig)

    # --- Plot 3: Individual high-quality UMAP for MP assignment ---
    if 'MP_assignment' in adata_umap.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata_umap, color='MP_assignment',
                   palette='tab20', size=40, alpha=0.8,
                   title='Meta-Program Assignment Distribution',
                   frameon=False, ax=ax, show=False,
                   legend_loc='right margin')
        plt.savefig(f"{save_path}/mp_assignment_umap_high_quality.pdf", dpi=300, bbox_inches='tight')
        plt.show()
        plt.close(fig)

    # --- Plot 4: Density plots for each MP ---
    if 'MP_assignment' in adata_umap.obs.columns:
        mp_categories_to_plot = [mp for mp in adata_umap.obs['MP_assignment'].unique() if mp not in ['Unresolved', '', None]]
        n_plots = len(mp_categories_to_plot)

        if n_plots > 0:
            # CORRECTION: The original subplot calculation was incorrect for odd numbers. This is the fix.
            ncols = 3
            nrows = int(np.ceil(n_plots / ncols))
            
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
            axes = axes.flatten()

            for i, mp in enumerate(mp_categories_to_plot):
                adata_umap.obs[f'is_{mp}'] = (adata_umap.obs['MP_assignment'] == mp).astype('category')
                sc.pl.umap(adata_umap, color=f'is_{mp}',
                           size=30, alpha=0.8, cmap='Reds',
                           title=f'{mp} Distribution',
                           frameon=False, ax=axes[i], show=False, legend_loc=None)

            for j in range(i + 1, len(axes)):
                axes[j].axis('off')
            
            plt.tight_layout()
            plt.savefig(f"{save_path}/mp_individual_density_umap.pdf", dpi=300, bbox_inches='tight')
            plt.show()
            plt.close(fig)

    print("  ✅ UMAP visualizations saved!")

    # --- Save UMAP coordinates ---
    print("  💾 Saving UMAP coordinates...")
    umap_coords = pd.DataFrame({
        'cell_barcode': adata_umap.obs.index,
        'UMAP_1': adata_umap.obsm['X_umap'][:, 0],
        'UMAP_2': adata_umap.obsm['X_umap'][:, 1],
        'MP_assignment': adata_umap.obs.get('MP_assignment', 'Unknown'),
        'patient': adata_umap.obs.get('patient', 'Unknown')
    })
    for col in mp_score_cols:
        if col in adata_umap.obs.columns:
            umap_coords[col] = adata_umap.obs[col].values
    
    umap_coords_path = f"{save_path}/umap_coordinates_with_mp_data.csv"
    umap_coords.to_csv(umap_coords_path, index=False)
    print(f"  💾 UMAP coordinates saved to: {umap_coords_path}")

    return adata_umap

def run_mp_comparison_analysis(adata, comparison_type, group_col, group1, group2,
                             save_path, description, paired_analysis=False, patient_col=None):
    """Generic function to run MP comparison analysis (unchanged)."""
    
    print(f"\n{description}")
    
    patient_mp_data = []
    
    for patient_id in adata.obs['patient'].unique():
        patient_data = adata[adata.obs['patient'] == patient_id]
        
        if paired_analysis:
            for stage in patient_data.obs[group_col].unique():
                stage_data = patient_data[patient_data.obs[group_col] == stage]
                if len(stage_data) < 5: continue
                
                mp_counts = stage_data.obs['MP_assignment'].value_counts()
                total_cells = len(stage_data)
                
                for mp in adata.obs['MP_assignment'].unique():
                    if mp != 'Unresolved':
                        fraction = mp_counts.get(mp, 0) / total_cells
                        patient_mp_data.append({'patient': patient_id, 'mp': mp, 'fraction': fraction, 'group': stage, 'total_cells': total_cells})
        else:
            if len(patient_data) < 10: continue
            
            group_value = patient_data.obs[group_col].iloc[0]
            mp_counts = patient_data.obs['MP_assignment'].value_counts()
            total_cells = len(patient_data)
            
            for mp in adata.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    fraction = mp_counts.get(mp, 0) / total_cells
                    patient_mp_data.append({'patient': patient_id, 'mp': mp, 'fraction': fraction, 'group': group_value, 'total_cells': total_cells})

    if not patient_mp_data:
        print("No valid patient data for analysis")
        return None
    
    patient_df = pd.DataFrame(patient_mp_data)
    mp_results = []
    
    for mp in patient_df['mp'].unique():
        mp_data = patient_df[patient_df['mp'] == mp]
        group1_fractions = mp_data[mp_data['group'] == group1]['fraction']
        group2_fractions = mp_data[mp_data['group'] == group2]['fraction']
        
        if len(group1_fractions) >= 3 and len(group2_fractions) >= 3:
            if paired_analysis:
                group1_patients = set(mp_data[mp_data['group'] == group1]['patient'])
                group2_patients = set(mp_data[mp_data['group'] == group2]['patient'])
                common_patients = group1_patients & group2_patients
                
                if len(common_patients) >= 3:
                    paired_group1 = [mp_data[(mp_data['patient'] == p) & (mp_data['group'] == group1)]['fraction'].iloc[0] for p in common_patients]
                    paired_group2 = [mp_data[(mp_data['patient'] == p) & (mp_data['group'] == group2)]['fraction'].iloc[0] for p in common_patients]
                    
                    t_stat, p_value = stats.ttest_rel(paired_group1, paired_group2)
                    differences = np.array(paired_group1) - np.array(paired_group2)
                    cohens_d = np.mean(differences) / np.std(differences, ddof=1) if np.std(differences, ddof=1) > 0 else 0
                    
                    mp_results.append({'mp': mp, 'group1_mean': np.mean(paired_group1), 'group2_mean': np.mean(paired_group2), 'n_paired_patients': len(common_patients), 't_statistic': t_stat, 'p_value': p_value, 'cohens_d': cohens_d, 'test_type': 'paired_t_test'})
            else:
                t_stat, p_value = stats.ttest_ind(group1_fractions, group2_fractions)
                pooled_std = np.sqrt(((len(group1_fractions) - 1) * np.var(group1_fractions, ddof=1) + (len(group2_fractions) - 1) * np.var(group2_fractions, ddof=1)) / (len(group1_fractions) + len(group2_fractions) - 2))
                cohens_d = (np.mean(group1_fractions) - np.mean(group2_fractions)) / pooled_std if pooled_std > 0 else 0
                
                mp_results.append({'mp': mp, 'group1_mean': np.mean(group1_fractions), 'group2_mean': np.mean(group2_fractions), 'n_group1_patients': len(group1_fractions), 'n_group2_patients': len(group2_fractions), 't_statistic': t_stat, 'p_value': p_value, 'cohens_d': cohens_d, 'test_type': 'unpaired_t_test'})

    if mp_results:
        mp_results_df = pd.DataFrame(mp_results)
        _, p_corrected, _, _ = multipletests(mp_results_df['p_value'], method='fdr_bh')
        mp_results_df['p_corrected'] = p_corrected
        mp_results_df['comparison_type'] = comparison_type
        mp_results_df['group1'] = group1
        mp_results_df['group2'] = group2
        
        filename = f"mp_{comparison_type}_resistance_analysis.csv"
        mp_results_df.to_csv(f"{save_path}/{filename}", index=False)
        return {'results_df': mp_results_df}
    
    return None

def calculate_treatment_induced_changes(adata, base_mask):
    """Calculate treatment-induced MP changes for each patient"""
    
    # Filter data
    analysis_data = adata[base_mask].copy()
    
    # Find patients with both Pre and Post samples
    patient_stages = analysis_data.obs.groupby('patient')['Treatment_Stage'].unique()
    paired_patients = []
    
    for patient_id, stages in patient_stages.items():
        if 'Pre' in stages and 'Post' in stages:
            paired_patients.append(patient_id)
    
    if len(paired_patients) == 0:
        return None
    
    treatment_changes = []
    
    for patient_id in paired_patients:
        patient_data = analysis_data[analysis_data.obs['patient'] == patient_id]
        response = patient_data.obs['Response'].iloc[0]
        
        pre_data = patient_data[patient_data.obs['Treatment_Stage'] == 'Pre']
        post_data = patient_data[patient_data.obs['Treatment_Stage'] == 'Post']
        
        if len(pre_data) >= 5 and len(post_data) >= 5:
            # Calculate MP fractions for Pre and Post
            for mp in analysis_data.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    pre_fraction = (pre_data.obs['MP_assignment'] == mp).mean()
                    post_fraction = (post_data.obs['MP_assignment'] == mp).mean()
                    
                    change = post_fraction - pre_fraction
                    
                    treatment_changes.append({
                        'patient': patient_id,
                        'mp': mp,
                        'pre_fraction': pre_fraction,
                        'post_fraction': post_fraction,
                        'change': change,
                        'response': response
                    })
    
    return pd.DataFrame(treatment_changes) if treatment_changes else None

def compare_treatment_changes(treatment_changes_df, save_path):
    """Compare treatment-induced changes between non-pCR vs pCR patients"""
    
    change_results = []
    
    for mp in treatment_changes_df['mp'].unique():
        mp_data = treatment_changes_df[treatment_changes_df['mp'] == mp]
        
        nonpcr_changes = mp_data[mp_data['response'] == 'non_pCR']['change']
        pcr_changes = mp_data[mp_data['response'] == 'pCR']['change']
        
        if len(nonpcr_changes) >= 3 and len(pcr_changes) >= 3:
            t_stat, p_value = stats.ttest_ind(nonpcr_changes, pcr_changes)
            
            # Calculate effect size
            pooled_std = np.sqrt(((len(nonpcr_changes) - 1) * np.var(nonpcr_changes, ddof=1) + 
                                (len(pcr_changes) - 1) * np.var(pcr_changes, ddof=1)) / 
                               (len(nonpcr_changes) + len(pcr_changes) - 2))
            
            cohens_d = (np.mean(nonpcr_changes) - np.mean(pcr_changes)) / pooled_std if pooled_std > 0 else 0
            
            change_results.append({
                'mp': mp,
                'nonpcr_change_mean': np.mean(nonpcr_changes),
                'pcr_change_mean': np.mean(pcr_changes),
                'nonpcr_change_std': np.std(nonpcr_changes),
                'pcr_change_std': np.std(pcr_changes),
                'n_nonpcr_patients': len(nonpcr_changes),
                'n_pcr_patients': len(pcr_changes),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d
            })
            
            print(f"\n{mp} (Treatment-induced change):")
            print(f"  Non-pCR change: {np.mean(nonpcr_changes):.3f} ± {np.std(nonpcr_changes):.3f}")
            print(f"  pCR change: {np.mean(pcr_changes):.3f} ± {np.std(pcr_changes):.3f}")
            print(f"  Difference: {np.mean(nonpcr_changes) - np.mean(pcr_changes):.3f}")
            print(f"  t-test: t={t_stat:.3f}, p={p_value:.3e}")
            print(f"  Cohen's d: {cohens_d:.3f}")
    
    if change_results:
        change_results_df = pd.DataFrame(change_results)
        
        # Multiple testing correction
        _, p_corrected, _, _ = multipletests(change_results_df['p_value'], method='fdr_bh')
        change_results_df['p_corrected'] = p_corrected
        
        # Save results
        change_results_df.to_csv(f"{save_path}/mp_therapy_specific_resistance.csv", index=False)
        
        # Identify significant therapy-resistant MPs
        therapy_resistant_mps = change_results_df[
            (change_results_df['p_corrected'] < 0.05) & 
            (change_results_df['cohens_d'] > 0.5) &  # More upregulated in non-pCR
            (change_results_df['nonpcr_change_mean'] > 0)  # Actually increased post-treatment
        ]
        
        if len(therapy_resistant_mps) > 0:
            print(f"\n🔥 THERAPY-RESISTANT META-PROGRAMS:")
            for _, row in therapy_resistant_mps.iterrows():
                print(f"  {row['mp']}: Upregulated in non-pCR after treatment")
                print(f"    Non-pCR change: +{row['nonpcr_change_mean']:.3f}")
                print(f"    pCR change: {row['pcr_change_mean']:.3f}")
                print(f"    Corrected p-value: {row['p_corrected']:.3e}")
        
        return {
            'results_df': change_results_df,
            'therapy_resistant_mps': therapy_resistant_mps,
            'comparison_type': 'therapy_specific'
        }
    
    return None

def create_resistance_summary(intrinsic_results, acquired_results, therapy_specific_results, save_path):
    """Create comprehensive resistance analysis summary"""
    
    print("\n📊 Creating Resistance Analysis Summary")
    
    # Combine all significant results
    all_significant_mps = []
    
    if intrinsic_results and 'significant_mps' in intrinsic_results:
        for _, row in intrinsic_results['significant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Intrinsic',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Enriched in pre-treatment non-pCR"
            })
    
    if acquired_results and 'significant_mps' in acquired_results:
        for _, row in acquired_results['significant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Acquired',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Upregulated post-treatment in non-pCR"
            })
    
    if therapy_specific_results and 'therapy_resistant_mps' in therapy_specific_results:
        for _, row in therapy_specific_results['therapy_resistant_mps'].iterrows():
            all_significant_mps.append({
                'mp': row['mp'],
                'resistance_type': 'Therapy-Specific',
                'effect_size': row['cohens_d'],
                'p_corrected': row['p_corrected'],
                'description': f"Differentially upregulated in non-pCR vs pCR"
            })
    
    if all_significant_mps:
        summary_df = pd.DataFrame(all_significant_mps)
        summary_df.to_csv(f"{save_path}/comprehensive_resistance_summary.csv", index=False)
        
        print(f"💡 Summary of Resistant Meta-Programs:")
        for resistance_type in summary_df['resistance_type'].unique():
            type_mps = summary_df[summary_df['resistance_type'] == resistance_type]
            print(f"  {resistance_type}: {len(type_mps)} MP(s)")
            for _, row in type_mps.iterrows():
                print(f"    - {row['mp']} (d={row['effect_size']:.3f}, p={row['p_corrected']:.3e})")
    
    print(f"\n✅ Comprehensive resistance analysis completed!")

def comprehensive_msi_vs_mss_resistance_analysis(adata, save_path):
    """
    Comprehensive resistance analysis comparing MSI vs MSS patterns
    
    Analysis strategy:
    1. MSI resistance patterns (previous analysis)
    2. MSS resistance patterns (new analysis)
    3. MSI vs MSS direct comparison
    4. Universal vs specific resistance mechanisms
    """
    
    print("\n🔬 Comprehensive MSI vs MSS Resistance Analysis")
    print("="*70)
    
    # Data overview across microsatellite status
    print("📊 Complete Data Overview:")
    total_treated_mask = adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx'])
    total_treated_data = adata[total_treated_mask].copy()
    
    print(f"Total treated cells: {total_treated_data.n_obs:,}")
    
    # Cross-tabulation of key variables
    print("\nMicrosatellite Status × Treatment Stage × Response:")
    detailed_crosstab = pd.crosstab(
        [total_treated_data.obs['Microsatellite_Status'], total_treated_data.obs['Treatment_Stage']], 
        total_treated_data.obs['Response'], 
        margins=True
    )
    print(detailed_crosstab)
    
    # Check MSS response rates
    mss_data = total_treated_data[total_treated_data.obs['Microsatellite_Status'] == 'MSS']
    if len(mss_data) > 0:
        print(f"\nMSS Response Analysis:")
        mss_response_counts = mss_data.obs['Response'].value_counts()
        print(mss_response_counts)
        
        if 'pCR' in mss_response_counts and 'non_pCR' in mss_response_counts:
            pCR_rate_mss = mss_response_counts['pCR'] / mss_response_counts.sum() * 100
            print(f"MSS pCR rate: {pCR_rate_mss:.1f}%")
        else:
            print("⚠️ MSS response data incomplete")
    
    # Analysis 1: MSI Resistance (Reference - already done)
    print("\n🎯 Analysis 1: MSI Resistance Patterns (Reference)")
    print("-" * 60)
    msi_results = analyze_microsatellite_specific_resistance(
        adata, 'MSI', save_path, "MSI Tumors"
    )
    
    # Analysis 2: MSS Resistance  
    print("\n🔍 Analysis 2: MSS Resistance Patterns")
    print("-" * 60)
    mss_results = analyze_microsatellite_specific_resistance(
        adata, 'MSS', save_path, "MSS Tumors"
    )
    
    # Analysis 3: Direct MSI vs MSS Comparison
    print("\n⚡ Analysis 3: Direct MSI vs MSS Comparison")
    print("-" * 60)
    direct_comparison_results = analyze_msi_vs_mss_programs(adata, save_path)
    
    # Analysis 4: Universal vs Specific Mechanisms
    print("\n🌐 Analysis 4: Universal vs Microsatellite-Specific Mechanisms")
    print("-" * 60)
    universal_vs_specific = identify_universal_vs_specific_resistance(
        msi_results, mss_results, save_path
    )
    
    # Analysis 5: Enhanced Clinical Insights
    print("\n🏥 Analysis 5: Enhanced Clinical Insights")
    print("-" * 60)
    clinical_insights = generate_clinical_insights(
        msi_results, mss_results, direct_comparison_results, save_path
    )
    
    return {
        'msi_resistance': msi_results,
        'mss_resistance': mss_results, 
        'direct_comparison': direct_comparison_results,
        'universal_vs_specific': universal_vs_specific,
        'clinical_insights': clinical_insights
    }

def analyze_microsatellite_specific_resistance(adata, msi_status, save_path, description):
    """Analyze resistance patterns within specific microsatellite status"""
    
    print(f"\n{description} Analysis")
    
    # Filter for specific microsatellite status
    ms_mask = (adata.obs['Microsatellite_Status'] == msi_status) & \
              (adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx']))
    
    ms_data = adata[ms_mask].copy()
    
    print(f"{msi_status} treated cells: {ms_data.n_obs:,}")
    
    if len(ms_data) == 0:
        print(f"⚠️ No {msi_status} treated cells found")
        return None
    
    # Check response distribution
    response_dist = ms_data.obs['Response'].value_counts()
    print(f"{msi_status} Response distribution:")
    print(response_dist)
    
    # Check if we have sufficient data for analysis
    if 'pCR' not in response_dist or 'non_pCR' not in response_dist:
        print(f"⚠️ Insufficient response data for {msi_status} analysis")
        return None
    
    if response_dist['pCR'] < 100 or response_dist['non_pCR'] < 100:
        print(f"⚠️ Low sample sizes for {msi_status}: pCR={response_dist.get('pCR', 0)}, non_pCR={response_dist.get('non_pCR', 0)}")
    
    # Run the three types of resistance analysis
    results = {}
    
    # 1. Intrinsic resistance (Pre-treatment)
    intrinsic_results = analyze_intrinsic_resistance_specific(ms_data, msi_status, save_path)
    results['intrinsic'] = intrinsic_results
    
    # 2. Acquired resistance (Post vs Pre in non-pCR) 
    acquired_results = analyze_acquired_resistance_specific(ms_data, msi_status, save_path)
    results['acquired'] = acquired_results
    
    # 3. Therapy-specific resistance
    therapy_results = analyze_therapy_specific_resistance_specific(ms_data, msi_status, save_path)
    results['therapy_specific'] = therapy_results
    
    return results

def analyze_intrinsic_resistance_specific(ms_data, msi_status, save_path):
    """Intrinsic resistance analysis for specific microsatellite status"""
    
    # Filter for pre-treatment
    pre_mask = ms_data.obs['Treatment_Stage'] == 'Pre'
    pre_data = ms_data[pre_mask].copy()
    
    print(f"\n{msi_status} Pre-treatment analysis:")
    print(f"  Cells: {pre_data.n_obs:,}")
    print(f"  Response distribution: {pre_data.obs['Response'].value_counts().to_dict()}")
    
    if len(pre_data) < 100:
        print(f"  ⚠️ Insufficient pre-treatment {msi_status} data")
        return None
    
    return run_mp_comparison_analysis(
        pre_data,
        comparison_type=f'{msi_status}_intrinsic',
        group_col='Response',
        group1='non_pCR',
        group2='pCR', 
        save_path=save_path,
        description=f"{msi_status} Pre-treatment non-pCR vs pCR",
        paired_analysis=False
    )

def analyze_acquired_resistance_specific(ms_data, msi_status, save_path):
    """Acquired resistance analysis for specific microsatellite status"""
    
    # Filter for non-pCR with paired samples
    nonpcr_data = ms_data[ms_data.obs['Response'] == 'non_pCR'].copy()
    
    # Find patients with both Pre and Post
    patient_stages = nonpcr_data.obs.groupby('patient')['Treatment_Stage'].unique()
    paired_patients = [pid for pid, stages in patient_stages.items() 
                      if 'Pre' in stages and 'Post' in stages]
    
    print(f"\n{msi_status} Acquired resistance analysis:")
    print(f"  Non-pCR patients with paired samples: {len(paired_patients)}")
    
    if len(paired_patients) < 3:
        print(f"  ⚠️ Insufficient paired {msi_status} data")
        return None
    
    paired_data = nonpcr_data[nonpcr_data.obs['patient'].isin(paired_patients)].copy()
    
    return run_mp_comparison_analysis(
        paired_data,
        comparison_type=f'{msi_status}_acquired',
        group_col='Treatment_Stage',
        group1='Post',
        group2='Pre',
        save_path=save_path,
        description=f"{msi_status} Post vs Pre in non-pCR",
        paired_analysis=True,
        patient_col='patient'
    )

def analyze_therapy_specific_resistance_specific(ms_data, msi_status, save_path):
    """Therapy-specific resistance for specific microsatellite status"""
    
    # Calculate treatment-induced changes
    ms_changes = calculate_treatment_induced_changes_ms(ms_data, msi_status)
    
    if ms_changes is None:
        print(f"  ⚠️ No treatment change data for {msi_status}")
        return None
    
    print(f"\n{msi_status} Therapy-specific analysis:")
    print(f"  Patients with change data: {len(ms_changes)}")
    
    return compare_treatment_changes_ms(ms_changes, msi_status, save_path)

def analyze_msi_vs_mss_programs(adata, save_path):
    """Direct comparison of MP patterns between MSI and MSS tumors"""
    
    print("\nDirect MSI vs MSS Meta-Program Comparison")
    
    # Filter for treated samples with clear microsatellite status
    treated_mask = adata.obs['Treatment_Strategy'].isin(['Anti-PD1', 'Anti-PD1 plus Celecoxib', 'Anti-PD1 plus CapeOx'])
    ms_clear_mask = adata.obs['Microsatellite_Status'].isin(['MSI', 'MSS'])
    
    comparison_data = adata[treated_mask & ms_clear_mask].copy()
    
    print(f"Comparison cohort: {comparison_data.n_obs:,} cells")
    print("Microsatellite distribution:")
    print(comparison_data.obs['Microsatellite_Status'].value_counts())
    
    # Overall MP distribution comparison
    overall_comparison = run_mp_comparison_analysis(
        comparison_data,
        comparison_type='msi_vs_mss_overall',
        group_col='Microsatellite_Status', 
        group1='MSS',
        group2='MSI',
        save_path=save_path,
        description="Overall MSS vs MSI Meta-Program Distribution",
        paired_analysis=False
    )
    
    # Pre-treatment specific comparison
    pre_comparison_data = comparison_data[comparison_data.obs['Treatment_Stage'] == 'Pre'].copy()
    
    if len(pre_comparison_data) > 100:
        pre_comparison = run_mp_comparison_analysis(
            pre_comparison_data,
            comparison_type='msi_vs_mss_pre',
            group_col='Microsatellite_Status',
            group1='MSS', 
            group2='MSI',
            save_path=save_path,
            description="Pre-treatment MSS vs MSI Meta-Program Patterns",
            paired_analysis=False
        )
    else:
        pre_comparison = None
    
    # Response-stratified comparison (within non-pCR and pCR separately)
    response_stratified = {}
    
    for response in ['pCR', 'non_pCR']:
        response_data = comparison_data[comparison_data.obs['Response'] == response].copy()
        
        if len(response_data) > 50:
            msi_count = (response_data.obs['Microsatellite_Status'] == 'MSI').sum()
            mss_count = (response_data.obs['Microsatellite_Status'] == 'MSS').sum()
            
            print(f"\n{response} stratified analysis: MSI={msi_count}, MSS={mss_count}")
            
            if msi_count >= 20 and mss_count >= 20:
                response_comp = run_mp_comparison_analysis(
                    response_data,
                    comparison_type=f'msi_vs_mss_{response}',
                    group_col='Microsatellite_Status',
                    group1='MSS',
                    group2='MSI', 
                    save_path=save_path,
                    description=f"MSS vs MSI in {response} patients",
                    paired_analysis=False
                )
                response_stratified[response] = response_comp
    
    return {
        'overall_comparison': overall_comparison,
        'pre_treatment_comparison': pre_comparison,
        'response_stratified': response_stratified
    }

def identify_universal_vs_specific_resistance(msi_results, mss_results, save_path):
    """Identify universal resistance mechanisms vs microsatellite-specific ones"""
    
    print("\nIdentifying Universal vs Specific Resistance Mechanisms")
    
    if not msi_results or not mss_results:
        print("⚠️ Insufficient data for universal vs specific analysis")
        return None
    
    # Extract significant MPs from each analysis
    msi_significant = extract_significant_mps(msi_results)
    mss_significant = extract_significant_mps(mss_results)
    
    # Find overlaps and differences
    universal_mps = []
    msi_specific_mps = []
    mss_specific_mps = []
    
    all_msi_mps = set(msi_significant.keys()) if msi_significant else set()
    all_mss_mps = set(mss_significant.keys()) if mss_significant else set()
    
    # Universal: significant in both MSI and MSS
    universal_mp_names = all_msi_mps & all_mss_mps
    
    # Specific: significant in only one
    msi_specific_mp_names = all_msi_mps - all_mss_mps
    mss_specific_mp_names = all_mss_mps - all_msi_mps
    
    print(f"\n🌐 Universal resistance MPs: {list(universal_mp_names)}")
    print(f"🎯 MSI-specific resistance MPs: {list(msi_specific_mp_names)}")
    print(f"🔍 MSS-specific resistance MPs: {list(mss_specific_mp_names)}")
    
    # Create summary
    universal_vs_specific_summary = {
        'universal_mps': universal_mp_names,
        'msi_specific_mps': msi_specific_mp_names,
        'mss_specific_mps': mss_specific_mp_names,
        'msi_results': msi_significant,
        'mss_results': mss_significant
    }
    
    # Save detailed results
    with open(f"{save_path}/universal_vs_specific_resistance.pkl", 'wb') as f:
        pickle.dump(universal_vs_specific_summary, f)
    
    return universal_vs_specific_summary

def extract_significant_mps(results_dict):
    """Extract significant MPs from resistance analysis results"""
    
    significant_mps = {}
    
    if not results_dict:
        return significant_mps
    
    for analysis_type, analysis_results in results_dict.items():
        if analysis_results and 'significant_mps' in analysis_results:
            for _, row in analysis_results['significant_mps'].iterrows():
                mp_name = row['mp']
                if mp_name not in significant_mps:
                    significant_mps[mp_name] = []
                
                significant_mps[mp_name].append({
                    'analysis_type': analysis_type,
                    'effect_size': row['cohens_d'],
                    'p_corrected': row['p_corrected']
                })
    
    return significant_mps

def generate_clinical_insights(msi_results, mss_results, direct_comparison, save_path):
    """Generate clinical insights from MSI vs MSS analysis"""
    
    print("\nGenerating Clinical Insights")
    
    insights = []
    
    # Response rate insights
    if msi_results and mss_results:
        insights.append("📊 Response Rate Patterns:")
        insights.append("   - Compare intrinsic resistance programs")
        insights.append("   - Identify why MSS typically doesn't respond")
        
    # Biomarker insights  
    if direct_comparison and 'overall_comparison' in direct_comparison:
        insights.append("🎯 Biomarker Implications:")
        insights.append("   - MSI-specific programs for current immunotherapy")
        insights.append("   - MSS-specific programs for future therapeutic targeting")
        
    # Therapeutic insights
    insights.append("💊 Therapeutic Insights:")
    insights.append("   - Universal targets: work across microsatellite status")
    insights.append("   - Specific targets: tailored to MSI or MSS tumors")
    
    # Save insights
    with open(f"{save_path}/clinical_insights.txt", 'w') as f:
        f.write("\n".join(insights))
    
    print("\n".join(insights))
    
    return insights

# Helper functions for MSS-specific analysis
def calculate_treatment_induced_changes_ms(ms_data, msi_status):
    """Calculate treatment changes for specific microsatellite status"""
    return calculate_treatment_induced_changes(ms_data, slice(None))  # Use existing function

def compare_treatment_changes_ms(changes_df, msi_status, save_path):
    """Compare treatment changes for specific microsatellite status"""
    # Modify the existing function to work with MS-specific data
    return compare_treatment_changes(changes_df, save_path)

# MSI-H CRC Specific Resistance Analysis Functions
def calculate_mp_fractions_per_sample(adata, mp_column='MP_assignment'):
    """
    Calculate MP fractions for each sample
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with cells
    mp_column : str
        Column name containing MP assignments
    
    Returns:
    --------
    pd.DataFrame
        Sample-level MP fractions with clinical metadata
    """
    
    print("📊 Calculating MP fractions per sample...")
    
    sample_mp_data = []
    
    for sample_id in adata.obs['sample_id'].unique():
        sample_cells = adata[adata.obs['sample_id'] == sample_id]
        
        if len(sample_cells) < 10:  # Minimum cells per sample
            continue
        
        # Get sample metadata (taking first cell's metadata as sample-level)
        sample_metadata = sample_cells.obs.iloc[0]
        
        # Calculate MP fractions
        mp_counts = sample_cells.obs[mp_column].value_counts()
        total_cells = len(sample_cells)
        
        # Get all unique MPs in the dataset
        all_mps = adata.obs[mp_column].unique()
        
        for mp in all_mps:
            if mp != 'Unresolved' and pd.notna(mp):
                fraction = mp_counts.get(mp, 0) / total_cells
                
                sample_mp_data.append({
                    'sample_id': sample_id,
                    'patient': sample_metadata.get('patient', sample_id),
                    'MP': mp,
                    'MP_fraction': fraction,
                    'total_cells': total_cells,
                    'mp_cell_count': mp_counts.get(mp, 0),
                    
                    # Clinical metadata
                    'Treatment_Strategy': sample_metadata.get('Treatment_Strategy', 'Unknown'),
                    'Microsatellite_Status': sample_metadata.get('Microsatellite_Status', 'Unknown'),
                    'Response': sample_metadata.get('Response', 'Unknown'),
                    'Treatment_Stage': sample_metadata.get('Treatment_Stage', 'Unknown'),
                    'Gender': sample_metadata.get('Gender', 'Unknown'),
                    'Age': sample_metadata.get('Age', np.nan),
                    'study': sample_metadata.get('study', 'Unknown'),
                    'Tissue_Type': sample_metadata.get('Tissue_Type', 'Unknown')
                })
    
    sample_df = pd.DataFrame(sample_mp_data)
    
    print(f"✅ Calculated MP fractions for {len(sample_df['sample_id'].unique())} samples")
    print(f"   Total MP-sample combinations: {len(sample_df)}")
    
    return sample_df

# Missing Key Figures Implementation
def create_nmf_robustness_plots(robust_results, save_path):
    """Figure 1: NMF Parameter Sensitivity Analysis"""
    
    print("📊 Creating NMF robustness plots...")
    
    # Extract robustness statistics
    all_programs = robust_results['all_programs']
    robust_programs = robust_results['robust_programs']
    parameters = robust_results['parameters']
    
    # Analyze program survival by K value
    k_survival = {}
    for k in parameters['K_range']:
        total_programs = len([p for p in all_programs if p['k_value'] == k])
        surviving_programs = len([p for p in robust_programs if p['k_value'] == k])
        k_survival[k] = {'total': total_programs, 'surviving': surviving_programs, 
                        'survival_rate': surviving_programs / total_programs if total_programs > 0 else 0}
    
    # Create subplot figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Program survival by K
    k_values = list(k_survival.keys())
    survival_rates = [k_survival[k]['survival_rate'] for k in k_values]
    total_programs = [k_survival[k]['total'] for k in k_values]
    
    axes[0,0].bar(k_values, survival_rates, alpha=0.7, color='steelblue')
    axes[0,0].set_xlabel('Number of NMF Components (K)')
    axes[0,0].set_ylabel('Program Survival Rate')
    axes[0,0].set_title('NMF Robustness: Program Survival by K Value')
    axes[0,0].set_ylim(0, 1)
    
    # Add text annotations
    for i, (k, rate) in enumerate(zip(k_values, survival_rates)):
        axes[0,0].text(k, rate + 0.02, f'{rate:.2f}', ha='center', fontweight='bold')
    
    # Plot 2: Total programs generated by K
    axes[0,1].bar(k_values, total_programs, alpha=0.7, color='orange')
    axes[0,1].set_xlabel('Number of NMF Components (K)')
    axes[0,1].set_ylabel('Total Programs Generated')
    axes[0,1].set_title('Program Generation by K Value')
    
    # Plot 3: Patient representation in robust programs
    patient_counts = {}
    for program in robust_programs:
        patient = program['patient']
        patient_counts[patient] = patient_counts.get(patient, 0) + 1
    
    patient_program_counts = list(patient_counts.values())
    axes[1,0].hist(patient_program_counts, bins=max(10, len(set(patient_program_counts))), 
                   alpha=0.7, color='green', edgecolor='black')
    axes[1,0].set_xlabel('Number of Robust Programs per Patient')
    axes[1,0].set_ylabel('Number of Patients')
    axes[1,0].set_title('Patient Contribution to Robust Programs')
    
    # Plot 4: Robustness criteria statistics
    criteria_names = ['Total Programs', 'Cross-K Stable', 'Cross-Patient', 'Final Robust']
    criteria_counts = [
        len(all_programs),
        len(all_programs),  # This would need to be tracked from the actual function
        len(all_programs),  # This would need to be tracked from the actual function  
        len(robust_programs)
    ]
    
    axes[1,1].bar(range(len(criteria_names)), criteria_counts, alpha=0.7, color='purple')
    axes[1,1].set_xticks(range(len(criteria_names)))
    axes[1,1].set_xticklabels(criteria_names, rotation=45, ha='right')
    axes[1,1].set_ylabel('Number of Programs')
    axes[1,1].set_title('Program Filtering Pipeline')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig1_nmf_robustness.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ NMF robustness plots saved")

def create_mp_assignment_quality_plots(malignant_cells, save_path):
    """Figure 2: MP Assignment Quality Assessment"""
    
    print("📊 Creating MP assignment quality plots...")
    
    # Get MP score columns
    mp_score_cols = [col for col in malignant_cells.obs.columns if 'MP' in col and 'score' in col and 'normalized' not in col]
    
    if len(mp_score_cols) == 0:
        print("⚠️ No MP scores found")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: MP score distributions
    score_data = malignant_cells.obs[mp_score_cols].values
    
    axes[0,0].boxplot(score_data, labels=[col.replace('_score', '') for col in mp_score_cols])
    axes[0,0].set_xlabel('Meta-Programs')
    axes[0,0].set_ylabel('MP Score')
    axes[0,0].set_title('MP Score Distributions')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Assignment confidence (max vs second-max score)
    max_scores = np.max(score_data, axis=1)
    second_max_scores = np.partition(score_data, -2, axis=1)[:, -2]
    confidence_scores = max_scores - second_max_scores
    
    axes[0,1].hist(confidence_scores, bins=50, alpha=0.7, color='orange', edgecolor='black')
    axes[0,1].axvline(x=np.median(confidence_scores), color='red', linestyle='--', 
                      label=f'Median: {np.median(confidence_scores):.3f}')
    axes[0,1].set_xlabel('Assignment Confidence (Max - Second Max Score)')
    axes[0,1].set_ylabel('Number of Cells')
    axes[0,1].set_title('MP Assignment Confidence Distribution')
    axes[0,1].legend()
    
    # Plot 3: MP assignment proportions
    if 'MP_assignment' in malignant_cells.obs.columns:
        assignment_counts = malignant_cells.obs['MP_assignment'].value_counts()
        
        # Create pie chart
        axes[1,0].pie(assignment_counts.values, labels=assignment_counts.index, autopct='%1.1f%%', startangle=90)
        axes[1,0].set_title('MP Assignment Distribution')
    
    # Plot 4: Score correlation matrix
    score_corr = malignant_cells.obs[mp_score_cols].corr()
    
    im = axes[1,1].imshow(score_corr, cmap='RdBu_r', vmin=-1, vmax=1)
    axes[1,1].set_xticks(range(len(mp_score_cols)))
    axes[1,1].set_yticks(range(len(mp_score_cols)))
    axes[1,1].set_xticklabels([col.replace('_score', '') for col in mp_score_cols], rotation=45)
    axes[1,1].set_yticklabels([col.replace('_score', '') for col in mp_score_cols])
    axes[1,1].set_title('MP Score Correlation Matrix')
    
    # Add correlation values as text
    for i in range(len(mp_score_cols)):
        for j in range(len(mp_score_cols)):
            text = axes[1,1].text(j, i, f'{score_corr.iloc[i, j]:.2f}', 
                                 ha="center", va="center", color="white" if abs(score_corr.iloc[i, j]) > 0.5 else "black")
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=axes[1,1], shrink=0.8)
    cbar.set_label('Correlation Coefficient')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig2_mp_assignment_quality.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ MP assignment quality plots saved")

def create_clustering_dendrogram(robust_results, save_path):
    """Figure 3: Hierarchical Clustering Dendrogram"""
    
    print("📊 Creating clustering dendrogram...")
    
    if 'clustering_info' not in robust_results:
        print("⚠️ No clustering info found")
        return
    
    clustering_info = robust_results['clustering_info']
    linkage_matrix = clustering_info['linkage_matrix']
    cluster_assignments = clustering_info['cluster_assignments']
    
    # Create dendrogram
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
    
    # Plot 1: Full dendrogram
    dendrogram(linkage_matrix, ax=ax1, color_threshold=0.7, above_threshold_color='gray')
    ax1.set_title('Hierarchical Clustering of Robust Programs')
    ax1.set_xlabel('Program Index')
    ax1.set_ylabel('Distance')
    
    # Plot 2: Cluster size distribution
    unique_clusters, cluster_sizes = np.unique(cluster_assignments, return_counts=True)
    
    bars = ax2.bar(range(len(unique_clusters)), cluster_sizes, alpha=0.7, color='steelblue')
    ax2.set_xlabel('Meta-Program ID')
    ax2.set_ylabel('Number of Programs')
    ax2.set_title('Programs per Meta-Program')
    ax2.set_xticks(range(len(unique_clusters)))
    ax2.set_xticklabels([f'MP{i}' for i in unique_clusters])
    
    # Add value labels on bars
    for bar, size in zip(bars, cluster_sizes):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1, f'{size}', 
                ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{save_path}/missing_fig3_clustering_dendrogram.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ Clustering dendrogram saved")


# Export function for R integration
def export_data_for_r_figures(malignant_cells, robust_results, save_path):
    """Export data needed for R figure generation"""
    
    print("📤 Exporting data for R figure generation...")
    
    # 1. Sample-level MP fractions for stacked barplot
    sample_mp_fractions = []
    
    for sample_id in malignant_cells.obs['sample_id'].unique():
        sample_cells = malignant_cells[malignant_cells.obs['sample_id'] == sample_id]
        sample_metadata = sample_cells.obs.iloc[0]
        
        if 'MP_assignment' in sample_cells.obs.columns:
            mp_counts = sample_cells.obs['MP_assignment'].value_counts()
            total_cells = len(sample_cells)
            
            for mp in malignant_cells.obs['MP_assignment'].unique():
                if mp != 'Unresolved':
                    fraction = mp_counts.get(mp, 0) / total_cells
                    
                    sample_mp_fractions.append({
                        'sample_id': sample_id,
                        'patient': sample_metadata.get('patient'),
                        'MP': mp,
                        'MP_fraction': fraction,
                        'MP_count': mp_counts.get(mp, 0),
                        'total_cells': total_cells,
                        'Microsatellite_Status': sample_metadata.get('Microsatellite_Status'),
                        'Response': sample_metadata.get('Response'),
                        'Treatment_Stage': sample_metadata.get('Treatment_Stage'),
                        'Treatment_Strategy': sample_metadata.get('Treatment_Strategy')
                    })
    
    sample_mp_df = pd.DataFrame(sample_mp_fractions)
    sample_mp_df.to_csv(f"{save_path}/r_sample_mp_fractions.csv", index=False)
    
    # 2. MP signature genes for heatmap
    if 'meta_programs' in robust_results:
        mp_signatures = {}
        all_signature_genes = set()
        
        for mp_id, mp_data in robust_results['meta_programs'].items():
            signature_genes = mp_data['signature_genes']
            mp_signatures[f'MP{mp_id}'] = signature_genes
            all_signature_genes.update(signature_genes)
        
        # Create binary signature matrix
        signature_matrix = pd.DataFrame(0, index=list(all_signature_genes), 
                                       columns=list(mp_signatures.keys()))
        
        for mp_name, genes in mp_signatures.items():
            for gene in genes:
                signature_matrix.loc[gene, mp_name] = 1
        
        signature_matrix.to_csv(f"{save_path}/r_mp_signature_matrix.csv")
        
        # Expression matrix for signature genes
        if all_signature_genes:
            available_genes = [g for g in all_signature_genes if g in malignant_cells.var_names]
            if available_genes:
                expr_subset = malignant_cells[:, available_genes]
                
                # Calculate mean expression per sample per gene
                sample_gene_expr = []
                for sample_id in malignant_cells.obs['sample_id'].unique():
                    sample_cells = expr_subset[expr_subset.obs['sample_id'] == sample_id]
                    mean_expr = sample_cells.X.mean(axis=0)
                    
                    if hasattr(mean_expr, 'A1'):  # Handle sparse matrices
                        mean_expr = mean_expr.A1
                    
                    for i, gene in enumerate(available_genes):
                        sample_gene_expr.append({
                            'sample_id': sample_id,
                            'Gene': gene,
                            'Mean_Expression': mean_expr[i]
                        })
                
                sample_expr_df = pd.DataFrame(sample_gene_expr)
                sample_expr_df.to_csv(f"{save_path}/r_sample_gene_expression.csv", index=False)
    
    # 3. MP similarity matrix
    if 'clustering_info' in robust_results:
        similarity_matrix = robust_results['clustering_info']['similarity_matrix']
        mp_names = [f'MP{i+1}' for i in range(len(similarity_matrix))]
        
        similarity_df = pd.DataFrame(similarity_matrix, index=mp_names, columns=mp_names)
        similarity_df.to_csv(f"{save_path}/r_mp_similarity_matrix.csv")
    
    print("✅ Data exported for R figure generation")
    print(f"📁 Files created:")
    print(f"   - r_sample_mp_fractions.csv (for stacked barplot)")
    print(f"   - r_mp_signature_matrix.csv (for signature heatmap)")
    print(f"   - r_sample_gene_expression.csv (for expression heatmap)")
    print(f"   - r_mp_similarity_matrix.csv (for MP correlation)")

    """Diagnose why clustering is problematic"""
    
    print("🔍 Diagnosing clustering issues...")
    
    robust_programs = robust_results.get('robust_programs', [])
    clustering_info = robust_results.get('clustering_info', {})
    
    if not robust_programs:
        print("❌ No robust programs found!")
        return
    
    print(f"📊 Analysis of {len(robust_programs)} robust programs:")
    
    # 1. Check program similarity distribution
    similarities = []
    for i in range(len(robust_programs)):
        for j in range(i+1, len(robust_programs)):
            genes1 = set(robust_programs[i]['genes'])
            genes2 = set(robust_programs[j]['genes'])
            jaccard = len(genes1 & genes2) / len(genes1 | genes2)
            similarities.append(jaccard)
    
    print(f"   Similarity stats:")
    print(f"   - Mean: {np.mean(similarities):.3f}")
    print(f"   - Median: {np.median(similarities):.3f}")
    print(f"   - Min: {np.min(similarities):.3f}")
    print(f"   - Max: {np.max(similarities):.3f}")
    print(f"   - 90th percentile: {np.percentile(similarities, 90):.3f}")
    
    # 2. Check if programs are too similar (causing massive clusters)
    high_similarity_count = sum(1 for s in similarities if s > 0.5)
    print(f"   - High similarity pairs (>0.5): {high_similarity_count}/{len(similarities)}")
    
    # 3. Patient distribution
    patient_counts = {}
    for program in robust_programs:
        patient = program.get('patient', 'Unknown')
        patient_counts[patient] = patient_counts.get(patient, 0) + 1
    
    print(f"   Programs per patient: {dict(patient_counts)}")
    
    # 4. K-factor distribution
    k_counts = {}
    for program in robust_programs:
        k = program.get('K', 'Unknown')
        k_counts[k] = k_counts.get(k, 0) + 1
    
    print(f"   Programs per K: {dict(k_counts)}")
    
    return {
        'similarities': similarities,
        'patient_counts': patient_counts,
        'k_counts': k_counts
    }

def calculate_mp_sample_distribution(adata, mp_names):
    """Calculate how many samples each MP appears in"""
    
    distribution = {}
    mp_names = ["MP"+ str(x) for x in mp_names]
    
    for mp_name in mp_names:
        if 'MP_assignment' in adata.obs.columns and 'sample_id' in adata.obs.columns:
            # Count samples where this MP is present
            mp_cells = adata.obs['MP_assignment'] == mp_name
            samples_with_mp = adata.obs.loc[mp_cells, 'sample_id'].nunique()
            distribution[mp_name] = samples_with_mp
        else:
            distribution[mp_name] = 0
    
    return distribution

def create_mp_clustering_heatmap(malignant_cells, robust_results, save_path, figsize=(20, 12)):
    """
    Create comprehensive heatmap showing MP clustering results and markers
    """
    
    print("🎨 Creating MP clustering and marker heatmap...")
    
    # Extract MP information
    meta_programs = robust_results['meta_programs']
    n_mps = len(meta_programs)
    
    # Prepare the figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 2, 1.2], width_ratios=[3, 0.1, 1], 
                         hspace=0.3, wspace=0.1)
    
    # =========================================================================
    # Panel 1: MP Score Heatmap (Top)
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    
    # Get MP scores and assignments
    mp_score_cols = [f'MP{i}_score' for i in range(1, n_mps + 1)]
    mp_scores = malignant_cells.obs[mp_score_cols].copy()
    mp_assignments = malignant_cells.obs['MP_assignment'].copy()
    
    # Sort cells by MP assignment
    sort_order = mp_assignments.argsort()
    mp_scores_sorted = mp_scores.iloc[sort_order]
    mp_assignments_sorted = mp_assignments.iloc[sort_order]
    
    # Create MP score heatmap
    mp_scores_normalized = mp_scores_sorted.T  # Transpose for MPs as rows
    
    im1 = ax1.imshow(mp_scores_normalized, aspect='auto', cmap='RdYlBu_r', 
                     vmin=mp_scores_normalized.min().min(), 
                     vmax=mp_scores_normalized.max().max())
    
    ax1.set_yticks(range(n_mps))
    ax1.set_yticklabels([f'MP{i}' for i in range(1, n_mps + 1)])
    ax1.set_xlabel('Cells (sorted by MP assignment)')
    ax1.set_title('Meta-Program Scores Across All Cells', fontsize=14, fontweight='bold')
    
    # Add MP assignment boundaries
    assignment_positions = []
    current_assignment = mp_assignments_sorted.iloc[0]
    current_pos = 0
    
    for i, assignment in enumerate(mp_assignments_sorted):
        if assignment != current_assignment:
            ax1.axvline(x=i-0.5, color='white', linewidth=2)
            assignment_positions.append((current_pos, i-1, current_assignment))
            current_assignment = assignment
            current_pos = i
    assignment_positions.append((current_pos, len(mp_assignments_sorted)-1, current_assignment))
    
    # Add assignment labels
    for start, end, assignment in assignment_positions:
        if end - start > 50:  # Only label if segment is large enough
            ax1.text((start + end) / 2, -0.8, assignment, ha='center', va='top', 
                    fontsize=10, rotation=45)
    
    # Add colorbar for MP scores
    cbar_ax1 = fig.add_subplot(gs[0, 1])
    cbar1 = plt.colorbar(im1, cax=cbar_ax1)
    cbar1.set_label('MP Score', rotation=270, labelpad=15)
    
    # =========================================================================
    # Panel 2: Signature Gene Expression Heatmap (Middle)
    # =========================================================================

    ax2 = fig.add_subplot(gs[1, 0])
    
    print("🔍 Selecting most discriminative signature genes...")
    
    # Enhanced discriminative gene selection
    def select_discriminative_genes(malignant_cells, meta_programs, mp_assignments, 
                                  max_genes_per_mp=8, min_specificity=0.3):
        """
        Select highly discriminative genes for each MP using specificity scoring
        """
        discriminative_genes = []
        gene_mp_labels = []
        gene_specificity_scores = []
        
        mp_ids = sorted(meta_programs.keys())
        
        for mp_id in mp_ids:
            mp_data = meta_programs[mp_id]
            candidate_genes = mp_data['signature_genes'][:20]  # Consider top 20
            
            # Filter available genes
            available_genes = [g for g in candidate_genes if g in malignant_cells.var_names]
            
            if not available_genes:
                continue
                
            # Calculate specificity for each gene
            gene_specificity = {}
            
            for gene in available_genes:
                # Get expression in current MP vs others
                current_mp_mask = mp_assignments == f'MP{mp_id}'
                other_mp_mask = (mp_assignments != f'MP{mp_id}') & (mp_assignments != 'Unresolved')
                
                if current_mp_mask.sum() == 0 or other_mp_mask.sum() == 0:
                    continue
                
                current_mp_expr = malignant_cells[current_mp_mask, gene].X
                other_mp_expr = malignant_cells[other_mp_mask, gene].X
                
                if hasattr(current_mp_expr, 'toarray'):
                    current_mp_expr = current_mp_expr.toarray().flatten()
                    other_mp_expr = other_mp_expr.toarray().flatten()
                
                # Calculate specificity score (mean_current / (mean_current + mean_others))
                mean_current = np.mean(current_mp_expr)
                mean_others = np.mean(other_mp_expr)
                
                if mean_current + mean_others > 0:
                    specificity = mean_current / (mean_current + mean_others)
                    
                    # Also consider fold change
                    fold_change = (mean_current + 1e-6) / (mean_others + 1e-6)
                    
                    # Combined score: specificity * log2(fold_change)
                    combined_score = specificity * np.log2(fold_change + 1)
                    gene_specificity[gene] = combined_score
            
            # Select top discriminative genes for this MP
            if gene_specificity:
                top_genes = sorted(gene_specificity.items(), key=lambda x: x[1], reverse=True)
                selected_genes = [gene for gene, score in top_genes 
                                if score > min_specificity][:max_genes_per_mp]
                
                discriminative_genes.extend(selected_genes)
                gene_mp_labels.extend([f'MP{mp_id}'] * len(selected_genes))
                gene_specificity_scores.extend([gene_specificity[gene] for gene in selected_genes])
                
                print(f"  MP{mp_id}: {len(selected_genes)} discriminative genes "
                      f"(avg specificity: {np.mean([gene_specificity[g] for g in selected_genes]):.3f})")
        
        return discriminative_genes, gene_mp_labels, gene_specificity_scores
    
    # Select discriminative genes
    all_signature_genes, mp_gene_labels, specificity_scores = select_discriminative_genes(
        malignant_cells, meta_programs, mp_assignments, max_genes_per_mp=8, min_specificity=0.2
    )
    
    if len(all_signature_genes) == 0:
        print("❌ No discriminative genes found!")
        return None
    
    print(f"📊 Selected {len(all_signature_genes)} discriminative genes total")
    
    # Get expression data for discriminative genes
    gene_expr = malignant_cells[:, all_signature_genes].X
    if hasattr(gene_expr, 'toarray'):
        gene_expr = gene_expr.toarray()
    
    # Sort by MP assignment (same order as Panel 1)
    gene_expr_sorted = gene_expr[sort_order, :]
    
    # Enhanced normalization: Z-score per gene across all cells
    gene_expr_zscore = np.zeros_like(gene_expr_sorted)
    for i in range(gene_expr_sorted.shape[1]):
        gene_values = gene_expr_sorted[:, i]
        if np.std(gene_values) > 0:
            gene_expr_zscore[:, i] = (gene_values - np.mean(gene_values)) / np.std(gene_values)
    
    # Transpose for plotting (genes as rows)
    gene_expr_plot = gene_expr_zscore.T
    
    # Create enhanced heatmap with better color scheme
    im2 = ax2.imshow(gene_expr_plot, aspect='auto', cmap='RdYlBu_r',
                     vmin=-2.5, vmax=2.5, interpolation='bilinear')
    
    # Enhanced gene labels with MP grouping
    ax2.set_yticks(range(len(all_signature_genes)))
    
    # Create enhanced gene labels with specificity info
    enhanced_labels = []
    for i, (gene, mp_label, spec_score) in enumerate(zip(all_signature_genes, mp_gene_labels, specificity_scores)):
        enhanced_labels.append(f"{gene} ({spec_score:.2f})")
    
    ax2.set_yticklabels(enhanced_labels, fontsize=8)
    ax2.set_xlabel('Cells (sorted by MP assignment)', fontsize=12)
    ax2.set_title(f'Discriminative Signature Genes (n={len(all_signature_genes)})', 
                  fontsize=14, fontweight='bold')
    
    # Add MP assignment boundaries with enhanced styling
    for i, (start, end, assignment) in enumerate(assignment_positions):
        # Alternating background colors for MP regions
        if i % 2 == 0:
            ax2.axvspan(start-0.5, end+0.5, alpha=0.1, color='gray', zorder=0)
        
        # Boundary lines
        if i > 0:  # Don't add line at the very beginning
            ax2.axvline(x=start-0.5, color='white', linewidth=2, alpha=0.8)
    
    # Enhanced MP color coding with separators
    mp_colors = plt.cm.Set2(np.linspace(0, 1, n_mps))  # Better color palette
    
    # Group genes by MP and add visual separators
    current_mp = None
    mp_boundaries = []
    
    for i, mp_label in enumerate(mp_gene_labels):
        mp_num = int(mp_label.replace('MP', '')) - 1
        color = mp_colors[mp_num % len(mp_colors)]
        
        # Add colored rectangle for MP identification
        rect = Rectangle((-len(all_signature_genes)*0.02, i-0.4), 
                        len(all_signature_genes)*0.015, 0.8, 
                        facecolor=color, alpha=0.9, clip_on=False, zorder=10)
        ax2.add_patch(rect)
        
        # Track MP boundaries for separators
        if current_mp != mp_label:
            if current_mp is not None:
                mp_boundaries.append(i-0.5)
            current_mp = mp_label
    
    # Add horizontal separators between MP gene groups
    for boundary in mp_boundaries:
        ax2.axhline(y=boundary, color='white', linewidth=2, alpha=0.8)
    
    # Enhanced colorbar
    cbar_ax2 = fig.add_subplot(gs[1, 1])
    cbar2 = plt.colorbar(im2, cax=cbar_ax2)
    cbar2.set_label('Z-score Expression', rotation=270, labelpad=15, fontsize=11)
    cbar2.set_ticks([-2, -1, 0, 1, 2])
    
    # Add MP legend on the left side
    legend_y_positions = []
    current_mp = None
    start_idx = 0
    
    for i, mp_label in enumerate(mp_gene_labels + [None]):  # Add None to trigger final group
        if mp_label != current_mp:
            if current_mp is not None:
                # Calculate center position for this MP group
                center_y = (start_idx + i - 1) / 2
                legend_y_positions.append((center_y, current_mp))
            current_mp = mp_label
            start_idx = i
    
    # Add MP labels on the left
    for y_pos, mp_label in legend_y_positions:
        mp_num = int(mp_label.replace('MP', '')) - 1
        color = mp_colors[mp_num % len(mp_colors)]
        
        # Add MP label
        ax2.text(-len(all_signature_genes)*0.04, y_pos, mp_label, 
                ha='right', va='center', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.7))
    
    # Add grid for better readability
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax2.set_axisbelow(True)
    
    print(f"✅ Enhanced discriminative heatmap created with {len(all_signature_genes)} genes")
    
    # =========================================================================
    # Panel 3: MP Assignment Distribution (Bottom)
    # =========================================================================
    ax3 = fig.add_subplot(gs[2, 0])
    
    # Count MP assignments
    assignment_counts = mp_assignments.value_counts()
    
    # Create bar plot
    bars = ax3.bar(range(len(assignment_counts)), assignment_counts.values, 
                   color=[mp_colors[i % n_mps] for i in range(len(assignment_counts))])
    
    ax3.set_xticks(range(len(assignment_counts)))
    ax3.set_xticklabels(assignment_counts.index, rotation=45, ha='right')
    ax3.set_ylabel('Number of Cells')
    ax3.set_title('Cell Distribution Across Meta-Programs', fontsize=14, fontweight='bold')
    
    # Add percentage labels on bars
    total_cells = len(mp_assignments)
    for i, (bar, count) in enumerate(zip(bars, assignment_counts.values)):
        percentage = count / total_cells * 100
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + total_cells*0.01,
                f'{percentage:.1f}%', ha='center', va='bottom', fontsize=10)
    
    # =========================================================================
    # Panel 4: MP Characteristics Summary (Right)
    # =========================================================================
    ax4 = fig.add_subplot(gs[:, 2])
    ax4.axis('off')
    
    # Create MP summary text
    summary_text = "Meta-Program Summary\n" + "="*25 + "\n\n"
    
    for mp_id in sorted(meta_programs.keys()):
        mp_data = meta_programs[mp_id]
        n_cells = sum(mp_assignments == f'MP{mp_id}')
        percentage = n_cells / total_cells * 100
        
        summary_text += f"MP{mp_id} ({percentage:.1f}%):\n"
        summary_text += f"  • {len(mp_data['signature_genes'])} signature genes\n"
        summary_text += f"  • {mp_data['n_programs']} constituent programs\n"
        summary_text += f"  • {len(mp_data['patients_represented'])} patients\n"
        
        # Show top 3 genes
        top_genes = mp_data['signature_genes'][:3]
        available_top_genes = [g for g in top_genes if g in malignant_cells.var_names]
        if available_top_genes:
            summary_text += f"  • Top genes: {', '.join(available_top_genes)}\n"
        summary_text += "\n"
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    # =========================================================================
    # Save and Display
    # =========================================================================
    plt.suptitle('Meta-Program Clustering Analysis Overview', fontsize=16, fontweight='bold', y=0.98)
    
    # Save figure
    output_path = f"{save_path}/mp_clustering_comprehensive_heatmap.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    
    print(f"✅ Comprehensive heatmap saved to: {output_path}")
    plt.show()
    
    return fig

def create_mp_signature_comparison_heatmap(malignant_cells, robust_results, save_path, figsize=(16, 10)):
    """
    Create focused heatmap comparing MP signature genes
    """
    
    print("🎨 Creating MP signature comparison heatmap...")
    
    meta_programs = robust_results['meta_programs']
    
    # Extract signature genes for each MP
    mp_signatures = {}
    for mp_id in sorted(meta_programs.keys()):
        mp_data = meta_programs[mp_id]
        # Get top 15 genes per MP for focused view
        signature_genes = mp_data['signature_genes'][:15]
        available_genes = [gene for gene in signature_genes if gene in malignant_cells.var_names]
        mp_signatures[f'MP{mp_id}'] = available_genes
    
    # Create combined signature gene list
    all_genes = []
    gene_mp_mapping = {}
    for mp, genes in mp_signatures.items():
        for gene in genes:
            if gene not in gene_mp_mapping:
                gene_mp_mapping[gene] = []
                all_genes.append(gene)
            gene_mp_mapping[gene].append(mp)
    
    print(f"Total unique signature genes: {len(all_genes)}")
    
    # Calculate average expression for each MP assignment
    mp_assignments = malignant_cells.obs['MP_assignment'].copy()
    unique_assignments = [f'MP{i}' for i in range(1, len(meta_programs) + 1)] + ['Unresolved']
    unique_assignments = [mp for mp in unique_assignments if mp in mp_assignments.values]
    
    # Create expression matrix: genes x MP assignments
    expr_matrix = np.zeros((len(all_genes), len(unique_assignments)))
    
    for j, assignment in enumerate(unique_assignments):
        assignment_mask = mp_assignments == assignment
        if assignment_mask.sum() > 0:
            assignment_cells = malignant_cells[assignment_mask]
            for i, gene in enumerate(all_genes):
                if gene in assignment_cells.var_names:
                    gene_expr = assignment_cells[:, gene].X
                    if hasattr(gene_expr, 'toarray'):
                        gene_expr = gene_expr.toarray().flatten()
                    expr_matrix[i, j] = np.mean(gene_expr)
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    # Normalize by row (z-score across MP assignments)
    expr_matrix_zscore = np.apply_along_axis(zscore, 1, expr_matrix)
    expr_matrix_zscore = np.nan_to_num(expr_matrix_zscore)  # Replace NaN with 0
    
    # Create heatmap
    im = ax.imshow(expr_matrix_zscore, aspect='auto', cmap='RdYlBu_r',
                   vmin=-2, vmax=2)
    
    # Set ticks and labels
    ax.set_xticks(range(len(unique_assignments)))
    ax.set_xticklabels(unique_assignments, rotation=45, ha='right')
    ax.set_yticks(range(len(all_genes)))
    ax.set_yticklabels(all_genes, fontsize=8)
    
    ax.set_xlabel('Meta-Program Assignment')
    ax.set_ylabel('Signature Genes')
    ax.set_title('Meta-Program Signature Gene Expression Comparison', 
                 fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Z-score Expression', rotation=270, labelpad=15)
    
    # Color-code genes by their primary MP
    mp_colors = plt.cm.Set3(np.linspace(0, 1, len(meta_programs)))
    for i, gene in enumerate(all_genes):
        primary_mp = gene_mp_mapping[gene][0]  # First MP this gene belongs to
        mp_num = int(primary_mp.replace('MP', '')) - 1
        color = mp_colors[mp_num % len(mp_colors)]
        ax.add_patch(Rectangle((-0.8, i-0.5), 0.3, 1, facecolor=color, alpha=0.8, clip_on=False))
    
    # Add MP legend
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=mp_colors[i], alpha=0.8, label=f'MP{i+1}') 
                      for i in range(len(meta_programs))]
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.15, 0.5))
    
    plt.tight_layout()
    
    # Save figure
    output_path = f"{save_path}/mp_signature_comparison_heatmap.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    
    print(f"✅ Signature comparison heatmap saved to: {output_path}")
    plt.show()
    
    return fig

def perform_pathway_enrichment_analysis(meta_programs, malignant_cells, 
                                      imc_aligned_pathways, background_genes=None):
    """
    Perform pathway enrichment analysis for each MP's signature genes
    """
    
    print("🧬 Performing pathway enrichment analysis for MPs...")
    
    if background_genes is None:
        # background_genes = list(malignant_cells.var_names) ## all genes in dataset

        # Use the union of all signature genes across meta-programs as background
        background_genes = []
        for mp_id in sorted(meta_programs.keys()):
            mp_data = meta_programs[mp_id]
            signature_genes = mp_data['signature_genes']
            background_genes.extend(signature_genes)
    
    enrichment_results = {}
    
    for mp_id in sorted(meta_programs.keys()):
        mp_data = meta_programs[mp_id]
        signature_genes = mp_data['signature_genes']
        
        # Filter for genes available in dataset
        available_signature_genes = [g for g in signature_genes if g in background_genes]
        
        if len(available_signature_genes) < 5:
            print(f"  ⚠️ MP{mp_id}: Too few genes ({len(available_signature_genes)}) for enrichment")
            continue
        
        mp_enrichments = {}
        
        for pathway_name, pathway_genes in imc_aligned_pathways.items():
            # Filter pathway genes available in dataset
            available_pathway_genes = [g for g in pathway_genes if g in background_genes]
            
            if len(available_pathway_genes) < 5:
                continue
            
            # Hypergeometric test
            # k: genes in both signature and pathway
            # M: total background genes
            # n: pathway genes in background
            # N: signature genes in background
            
            overlap_genes = set(available_signature_genes) & set(available_pathway_genes)
            k = len(overlap_genes)  # successes in sample
            M = len(background_genes)  # population size
            n = len(available_pathway_genes)  # successes in population
            N = len(available_signature_genes)  # sample size
            
            # P-value: probability of observing k or more overlaps by chance
            p_value = hypergeom.sf(k-1, M, n, N)
            
            # Enrichment score (observed/expected)
            expected = (N * n) / M
            enrichment_score = k / expected if expected > 0 else 0
            
            # Only keep significant enrichments
            if p_value < 0.05 and k >= 2:
                mp_enrichments[pathway_name] = {
                    'overlap_genes': list(overlap_genes),
                    'overlap_count': k,
                    'pathway_size': n,
                    'signature_size': N,
                    'p_value': p_value,
                    'enrichment_score': enrichment_score,
                    'adjusted_p_value': p_value * len(imc_aligned_pathways)  # Bonferroni
                }
        
        if mp_enrichments:
            enrichment_results[f'MP{mp_id}'] = mp_enrichments
            print(f"  ✅ MP{mp_id}: {len(mp_enrichments)} significantly enriched pathways")
        else:
            print(f"  ⚠️ MP{mp_id}: No significantly enriched pathways")
    
    return enrichment_results

def create_mp_pathway_comprehensive_figure(malignant_cells, robust_results, 
                                         imc_aligned_pathways, save_path, 
                                         figsize=(24, 16)):
    """
    Create comprehensive MP-pathway relationship figure with improved layout
    """
    
    print("🎨 Creating comprehensive MP-pathway relationship figure...")
    
    # Perform enrichment analysis
    enrichment_results = perform_pathway_enrichment_analysis(
        robust_results['meta_programs'], malignant_cells, imc_aligned_pathways
    )
    
    if not enrichment_results:
        print("❌ No enrichment results found!")
        return None
    
    # Prepare the figure with better layout
    fig = plt.figure(figsize=figsize)
    
    # Create a 4x4 grid layout for better organization
    gs = fig.add_gridspec(4, 4, 
                         height_ratios=[0.1, 2, 1.2, 0.5], 
                         width_ratios=[2.5, 1, 1.2, 0.3],
                         hspace=0.25, wspace=0.15,
                         left=0.06, right=0.94, top=0.92, bottom=0.08)
    
    # =========================================================================
    # Title Section
    # =========================================================================
    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis('off')
    title_ax.text(0.5, 0.5, 'Meta-Program Pathway Enrichment Analysis', 
                  ha='center', va='center', fontsize=20, fontweight='bold',
                  transform=title_ax.transAxes)
    
    # =========================================================================
    # Panel 1: Main Enrichment Heatmap (Top Left - Large)
    # =========================================================================
    ax1 = fig.add_subplot(gs[1, :2])
    
    # Create enrichment matrix
    all_pathways = set()
    for mp_enrichments in enrichment_results.values():
        all_pathways.update(mp_enrichments.keys())
    
    all_pathways = sorted(list(all_pathways))
    all_mps = sorted(enrichment_results.keys())
    
    # Create matrices for heatmap
    enrichment_matrix = np.zeros((len(all_pathways), len(all_mps)))
    pvalue_matrix = np.ones((len(all_pathways), len(all_mps)))
    
    for j, mp in enumerate(all_mps):
        for i, pathway in enumerate(all_pathways):
            if pathway in enrichment_results[mp]:
                enrichment_matrix[i, j] = enrichment_results[mp][pathway]['enrichment_score']
                pvalue_matrix[i, j] = enrichment_results[mp][pathway]['p_value']
    
    # Create the main heatmap
    # Use log2 transform for enrichment scores for better visualization
    enrichment_matrix_log = np.log2(enrichment_matrix + 1)
    
    im1 = ax1.imshow(enrichment_matrix_log, aspect='auto', cmap='Reds', 
                     vmin=0, vmax=np.percentile(enrichment_matrix_log[enrichment_matrix_log > 0], 95))
    
    # Add significance markers
    for i in range(len(all_pathways)):
        for j in range(len(all_mps)):
            if pvalue_matrix[i, j] < 0.001:
                marker = '***'
            elif pvalue_matrix[i, j] < 0.01:
                marker = '**'
            elif pvalue_matrix[i, j] < 0.05:
                marker = '*'
            else:
                marker = ''
            
            if marker and enrichment_matrix[i, j] > 1:
                ax1.text(j, i, marker, ha='center', va='center', 
                        fontsize=10, fontweight='bold', color='white')
                
    # Format axes
    ax1.set_xticks(range(len(all_mps)))
    ax1.set_xticklabels(all_mps, rotation=0, ha='center', fontsize=12, fontweight='bold')
    ax1.set_yticks(range(len(all_pathways)))
    
    # Clean pathway names for display
    clean_pathway_names = [p.replace('HALLMARK_', '').replace('_', ' ').title() 
                          for p in all_pathways]
    ax1.set_yticklabels(clean_pathway_names, fontsize=10)
    
    ax1.set_title('Pathway Enrichment Heatmap', 
                  fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel('Meta-Programs', fontsize=12, fontweight='bold')
    ax1.set_ylabel('HALLMARK Pathways', fontsize=12, fontweight='bold')
    
    # Add grid for better readability
    ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax1.set_axisbelow(True)
    
    # =========================================================================
    # Panel 2: Network View (Top Right)
    # =========================================================================
    ax2 = fig.add_subplot(gs[1, 2])
    
    # Create network graph
    G = nx.Graph()
    
    # Add MP nodes
    mp_colors = plt.cm.Set2(np.linspace(0, 1, len(all_mps)))
    pos = {}
    
    # Position MPs in a circle
    mp_angles = np.linspace(0, 2*np.pi, len(all_mps), endpoint=False)
    mp_radius = 1.5
    
    for i, mp in enumerate(all_mps):
        x = mp_radius * np.cos(mp_angles[i])
        y = mp_radius * np.sin(mp_angles[i])
        pos[mp] = (x, y)
        G.add_node(mp, node_type='MP', color=mp_colors[i])
    
    # Add pathway nodes and edges (only top enriched pathways for clarity)
    top_pathways = []
    for mp in all_mps:
        if mp in enrichment_results:
            mp_enrichments = enrichment_results[mp]
            top_mp_pathways = sorted(mp_enrichments.items(), 
                                   key=lambda x: x[1]['enrichment_score'], reverse=True)[:3]
            for pathway, data in top_mp_pathways:
                if data['p_value'] < 0.05 and data['enrichment_score'] > 1.5:
                    top_pathways.append(pathway)
    
    top_pathways = list(set(top_pathways))  # Remove duplicates
    pathway_angles = np.linspace(0, 2*np.pi, len(top_pathways), endpoint=False)
    pathway_radius = 0.8
    
    for i, pathway in enumerate(top_pathways):
        clean_name = pathway.replace('HALLMARK_', '').replace('_', '\n')[:15] + '...' if len(pathway) > 15 else pathway.replace('HALLMARK_', '').replace('_', '\n')
        x = pathway_radius * np.cos(pathway_angles[i])
        y = pathway_radius * np.sin(pathway_angles[i])
        pos[clean_name] = (x, y)
        G.add_node(clean_name, node_type='pathway', color='lightblue')
        
        # Add edges for significant enrichments
        for j, mp in enumerate(all_mps):
            if mp in enrichment_results and pathway in enrichment_results[mp]:
                enrichment_score = enrichment_results[mp][pathway]['enrichment_score']
                p_value = enrichment_results[mp][pathway]['p_value']
                
                if p_value < 0.05 and enrichment_score > 1.5:
                    G.add_edge(mp, clean_name, weight=enrichment_score)
    
    # Draw network
    mp_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'MP']
    pathway_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'pathway']
    
    # Draw MP nodes
    nx.draw_networkx_nodes(G, pos, nodelist=mp_nodes, 
                          node_color=[G.nodes[node]['color'] for node in mp_nodes],
                          node_size=1000, ax=ax2, alpha=0.8)
    
    # Draw pathway nodes
    nx.draw_networkx_nodes(G, pos, nodelist=pathway_nodes,
                          node_color='lightblue', node_size=400, ax=ax2, alpha=0.6)
    
    # Draw edges with thickness proportional to enrichment
    edges = G.edges()
    if edges:
        weights = [G[u][v]['weight'] for u, v in edges]
        nx.draw_networkx_edges(G, pos, edgelist=edges, 
                              width=[max(1, w/2) for w in weights], alpha=0.6, ax=ax2)
    
    # Add labels
    nx.draw_networkx_labels(G, pos, labels={node: node for node in mp_nodes}, 
                           font_size=9, font_weight='bold', ax=ax2)
    nx.draw_networkx_labels(G, pos, labels={node: node for node in pathway_nodes}, 
                           font_size=6, ax=ax2)
    
    ax2.set_title('MP-Pathway Network\n(Top Enrichments)', fontsize=12, fontweight='bold')
    ax2.axis('off')
    
    # =========================================================================
    # Panel 3: Colorbar (Top Far Right)
    # =========================================================================
    cbar_ax = fig.add_subplot(gs[1, 3])
    cbar1 = plt.colorbar(im1, cax=cbar_ax)
    cbar1.set_label('log₂(Enrichment Score + 1)', rotation=270, labelpad=20, fontsize=11)
    
    # =========================================================================
    # Panel 4: Top Pathways Summary (Middle Left)
    # =========================================================================
    ax3 = fig.add_subplot(gs[2, :2])
    ax3.axis('off')
    
    # Create a formatted table-like summary
    summary_data = []
    headers = ['Meta-Program', 'Top Pathway', 'Enrichment Score', 'P-value']
    
    for mp in all_mps:
        if mp in enrichment_results:
            mp_enrichments = enrichment_results[mp]
            
            # Get top pathway
            if mp_enrichments:
                top_pathway = max(mp_enrichments.items(), 
                                key=lambda x: x[1]['enrichment_score'])
                pathway_name, data = top_pathway
                clean_name = pathway_name.replace('HALLMARK_', '').replace('_', ' ').title()
                
                summary_data.append([
                    mp,
                    clean_name[:30] + '...' if len(clean_name) > 30 else clean_name,
                    f"{data['enrichment_score']:.2f}",
                    f"{data['p_value']:.2e}"
                ])
    
    # Create table
    table = ax3.table(cellText=summary_data, 
                     colLabels=headers,
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.15, 0.5, 0.15, 0.2])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style the table
    for i in range(len(headers)):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    for i in range(1, len(summary_data) + 1):
        for j in range(len(headers)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')
    
    ax3.set_title('Top Enriched Pathway per Meta-Program', 
                  fontsize=14, fontweight='bold', pad=20)
    
    # =========================================================================
    # Panel 5: IMC Alignment (Middle Right)
    # =========================================================================
    ax4 = fig.add_subplot(gs[2, 2:])
    ax4.axis('off')
    
    # Add pathway-IMC marker mapping
    imc_mapping_text = "Pathway → IMC Marker Alignment\n" + "="*35 + "\n\n"
    pathway_imc_mapping = {
        'HALLMARK_GLYCOLYSIS': 'GLUT1, HK2',
        'HALLMARK_HYPOXIA': 'CA-IX, VEGF', 
        'HALLMARK_FATTY_ACID_METABOLISM': 'FASN',
        'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': 'Vimentin, EpCAM',
        'HALLMARK_E2F_TARGETS': 'Ki67',
        'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'COX4',
        'HALLMARK_MTORC1_SIGNALING': 'mTOR, S6',
        'HALLMARK_INFLAMMATORY_RESPONSE': 'CD68, CD163'
    }
    
    for pathway, markers in pathway_imc_mapping.items():
        clean_pathway = pathway.replace('HALLMARK_', '').replace('_', ' ').title()
        imc_mapping_text += f"• {clean_pathway}:\n  {markers}\n\n"
    
    ax4.text(0.05, 0.95, imc_mapping_text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor='lightyellow', alpha=0.8))
    
    # =========================================================================
    # Panel 6: Legend and Statistics (Bottom)
    # =========================================================================
    ax5 = fig.add_subplot(gs[3, :])
    ax5.axis('off')
    
    # Create comprehensive legend and statistics
    stats_text = f"Analysis Summary: {len(all_mps)} Meta-Programs • {len(all_pathways)} Significantly Enriched Pathways • {len(imc_aligned_pathways)} Total Pathways Analyzed"
    
    legend_text = "Statistical Significance: *** p<0.001, ** p<0.01, * p<0.05 | Enrichment Score = Observed/Expected Gene Overlap"
    
    method_text = "Methods: Hypergeometric test for pathway enrichment | Network shows top 3 pathways per MP (p<0.05, ES>1.5)"
    
    # Add text elements with different styling
    ax5.text(0.5, 0.75, stats_text, transform=ax5.transAxes, fontsize=12,
             ha='center', va='center', fontweight='bold',
             bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.8))
    
    ax5.text(0.5, 0.45, legend_text, transform=ax5.transAxes, fontsize=10,
             ha='center', va='center',
             bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.8))
    
    ax5.text(0.5, 0.15, method_text, transform=ax5.transAxes, fontsize=9,
             ha='center', va='center', style='italic',
             bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8))
    
    # =========================================================================
    # Save and Display
    # =========================================================================
    
    # Save figure
    output_path = f"{save_path}/mp_pathway_enrichment_analysis.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf', 
                facecolor='white', edgecolor='none')
    
    print(f"✅ MP-pathway figure saved to: {output_path}")
    plt.show()
    
    return fig, enrichment_results

# Usage function
def create_mp_pathway_analysis(malignant_cells, robust_results, save_path):
    """
    Main function to create MP-pathway analysis
    """
    
    # Define IMC-aligned HALLMARK pathways (use your function from before)
    imc_aligned_pathways = get_imc_aligned_hallmark_pathways() # get_all_hallmark_pathways()
    
    # Create comprehensive figure
    fig, enrichment_results = create_mp_pathway_comprehensive_figure(
        malignant_cells, robust_results, imc_aligned_pathways, save_path
    )
    
    # Print summary
    print("\n📊 MP-Pathway Enrichment Summary:")
    print("="*50)
    
    for mp, pathways in enrichment_results.items():
        print(f"\n{mp} ({len(pathways)} enriched pathways):")
        top_pathways = sorted(pathways.items(), 
                            key=lambda x: x[1]['enrichment_score'], reverse=True)[:3]
        
        for pathway, data in top_pathways:
            clean_name = pathway.replace('HALLMARK_', '').replace('_', ' ')
            print(f"  • {clean_name}: ES={data['enrichment_score']:.2f}, p={data['p_value']:.2e}")
    
    return fig, enrichment_results

def calculate_mp_scores(bulk_exp: pd.DataFrame, mp_signatures: dict) -> pd.DataFrame:
    """
    Calculates single-sample Gene Set Enrichment Analysis (ssGSEA) scores for
    each meta-program in each sample of a bulk RNA-seq dataset.

    Args:
        bulk_exp (pd.DataFrame): DataFrame of bulk expression data, with genes
                                 as rows and samples as columns.
        mp_signatures (dict): Dictionary where keys are meta-program IDs and
                              values are lists of signature genes.

    Returns:
        pd.DataFrame: A DataFrame containing the ssGSEA enrichment scores, with
                      samples as rows and meta-programs as columns.
    """
    print("Calculating ssGSEA scores for each meta-program...")
    # Run ssGSEA using gseapy
    ssgsea_results = gp.ssgsea(
        data=bulk_exp,
        gene_sets=mp_signatures,
        scale=True,  # Normalizes scores
        verbose=False
    )
    
    # The result has MPs as rows and samples as columns, so we transpose it
    mp_scores = ssgsea_results.res2d
    mp_scores.columns = ['Sample_ID',"MP","ES","NES"]
    print("✅ ssGSEA scores calculated successfully.")
    return mp_scores

def find_optimal_cutpoint(
    scores: pd.Series,
    time_data: pd.Series,
    status_data: pd.Series,
    min_group_size: int = 10,
    percentile_range: tuple = (10, 90)
):
    """
    Find optimal cutpoint for survival analysis using log-rank test.
    
    Args:
        scores (pd.Series): Continuous variable scores
        time_data (pd.Series): Survival times
        status_data (pd.Series): Event status (1=event, 0=censored)
        min_group_size (int): Minimum size for each group
        percentile_range (tuple): Range of percentiles to test (min_percentile, max_percentile)
    
    Returns:
        dict: Dictionary containing optimal cutpoint and related statistics
    """
    # Remove any missing values
    valid_idx = ~(scores.isna() | time_data.isna() | status_data.isna())
    scores_clean = scores[valid_idx]
    time_clean = time_data[valid_idx]
    status_clean = status_data[valid_idx]
    
    if len(scores_clean) < 2 * min_group_size:
        raise ValueError(f"Not enough samples to perform optimal cutpoint analysis. "
                        f"Need at least {2 * min_group_size} samples, got {len(scores_clean)}")
    
    # Generate candidate cutpoints between specified percentiles
    min_percentile, max_percentile = percentile_range
    min_cutpoint = np.percentile(scores_clean, min_percentile)
    max_cutpoint = np.percentile(scores_clean, max_percentile)
    
    # Create candidate cutpoints
    sorted_scores = np.sort(scores_clean.unique())
    candidate_cutpoints = sorted_scores[
        (sorted_scores >= min_cutpoint) & (sorted_scores <= max_cutpoint)
    ]
    
    best_pvalue = 1.0
    best_cutpoint = None
    best_stats = None
    cutpoint_results = []
    
    for cutpoint in candidate_cutpoints:
        # Split into groups
        high_mask = scores_clean > cutpoint
        low_mask = scores_clean <= cutpoint
        
        # Check minimum group sizes
        if np.sum(high_mask) < min_group_size or np.sum(low_mask) < min_group_size:
            continue
        
        try:
            # Perform log-rank test
            results = logrank_test(
                durations_A=time_clean[high_mask],
                durations_B=time_clean[low_mask],
                event_observed_A=status_clean[high_mask],
                event_observed_B=status_clean[low_mask]
            )
            
            p_value = results.p_value
            
            # Store results
            cutpoint_results.append({
                'cutpoint': cutpoint,
                'p_value': p_value,
                'test_statistic': results.test_statistic,
                'high_group_size': np.sum(high_mask),
                'low_group_size': np.sum(low_mask)
            })
            
            # Update best cutpoint
            if p_value < best_pvalue:
                best_pvalue = p_value
                best_cutpoint = cutpoint
                best_stats = {
                    'test_statistic': results.test_statistic,
                    'high_group_size': np.sum(high_mask),
                    'low_group_size': np.sum(low_mask)
                }
                
        except Exception as e:
            # Skip this cutpoint if log-rank test fails
            continue
    
    if best_cutpoint is None:
        raise ValueError("Could not find any valid cutpoint. Try adjusting min_group_size or percentile_range.")
    
    # Apply multiple testing correction (Bonferroni)
    n_tests = len(cutpoint_results)
    corrected_pvalue = min(best_pvalue * n_tests, 1.0)
    
    return {
        'optimal_cutpoint': best_cutpoint,
        'p_value': best_pvalue,
        'corrected_p_value': corrected_pvalue,
        'n_tests': n_tests,
        'test_statistic': best_stats['test_statistic'],
        'high_group_size': best_stats['high_group_size'],
        'low_group_size': best_stats['low_group_size'],
        'all_results': pd.DataFrame(cutpoint_results)
    }

def perform_survival_analysis(
    mp_scores: pd.DataFrame,
    clinical_data: pd.DataFrame,
    mp_id: str,
    time_col: str = 'OS_time',
    status_col: str = 'OS_status',
    cutpoint_method: str = 'optimal',
    save_path: str = 'results/survival_analysis'
):
    """
    Performs Kaplan-Meier survival analysis for a single meta-program.

    Args:
        mp_scores (pd.DataFrame): DataFrame of MP enrichment scores (samples x MPs).
        clinical_data (pd.DataFrame): DataFrame with clinical data, including
                                      survival time and status. Must have a
                                      'Sample_ID' column to merge.
        mp_id (str): The identifier for the meta-program to analyze (must be a
                     column name in mp_scores).
        time_col (str): Column name for survival time in clinical_data.
        status_col (str): Column name for survival status (1=event, 0=censored).
    """

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Merge clinical data with the specific MP scores
    # Ensure the index of mp_scores (Sample_ID) can be merged with the clinical data
    mp_scores = mp_scores.loc[mp_scores["MP"] == mp_id,]

    clinical_data['Sample_ID'] = clinical_data['Sample_ID'].astype(str)
    mp_scores['Sample_ID'] = mp_scores['Sample_ID']

    # Merge the DataFrames on 'Sample_ID'
    merged_data = clinical_data.merge(mp_scores, on='Sample_ID', how='inner')

    # Stratify patients into 'High' and 'Low' groups based on the median score
    score = merged_data['NES']

    # Determine cutpoint based on method
    if cutpoint_method == 'median':
        cutpoint = score.median()
        cutpoint_info = {
            'cutpoint': cutpoint,
            'method': 'median',
            'high_group_size': sum(score > cutpoint),
            'low_group_size': sum(score <= cutpoint)
        }
        
    elif cutpoint_method == 'optimal':
        optimal_results = find_optimal_cutpoint(
            scores=score,
            time_data=merged_data[time_col],
            status_data=merged_data[status_col],
            min_group_size=5,
            percentile_range=(10, 90)
        )
        cutpoint = optimal_results['optimal_cutpoint']
        cutpoint_info = optimal_results.copy()
        cutpoint_info['method'] = 'optimal'

    # Stratify patients into 'High' and 'Low' groups
    merged_data['Group'] = ['High' if x > cutpoint else 'Low' for x in score]
    
    # Separate data for each group
    high_group = merged_data[merged_data['Group'] == 'High']
    low_group = merged_data[merged_data['Group'] == 'Low']

    # Perform log-rank test
    results = logrank_test(
        durations_A=high_group[time_col],
        durations_B=low_group[time_col],
        event_observed_A=high_group[status_col],
        event_observed_B=low_group[status_col]
    )
    p_value = results.p_value

    # Plot Kaplan-Meier curves
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(7, 5))

    kmf.fit(low_group[time_col], low_group[status_col], label=f'Low Expression (n={len(low_group)})')
    kmf.plot_survival_function(ax=ax)

    kmf.fit(high_group[time_col], high_group[status_col], label=f'High Expression (n={len(high_group)})')
    kmf.plot_survival_function(ax=ax)
    
    # Add plot details
    ax.set_title(f'Survival Analysis for Meta-Program {mp_id}')
    ax.set_xlabel('Time (e.g., Days)')
    ax.set_ylabel('Overall Survival Probability')
    ax.text(
        0.05, 0.05,
        f'Log-rank p-value: {p_value:.4f}',
        transform=ax.transAxes,
        fontsize=12,
        bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5)
    )
    plt.tight_layout()
    output_path = f"{save_path}/survival_analysis_{mp_id}.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.show()
    
    print(f"Analysis for MP {mp_id}: Log-rank p-value = {p_value:.4f}")
    return p_value