import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc

from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches


def explore_clinical_data(adata, variables=None, show_counts=True, show_percentages=True):
    """
    Explore clinical variables in AnnData object
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    variables : list, optional
        List of variables to explore. If None, explores all categorical variables
    show_counts : bool
        Whether to show raw counts
    show_percentages : bool
        Whether to show percentages
    
    Returns:
    --------
    dict : Dictionary containing summary statistics for each variable
    """
    
    if variables is None:
        # Auto-detect categorical and object columns
        variables = []
        for col in adata.obs.columns:
            if adata.obs[col].dtype in ['object', 'category'] or adata.obs[col].dtype.name == 'category':
                variables.append(col)
    
    results = {}
    
    print("="*80)
    print(f"CLINICAL DATA EXPLORATION - Total cells: {adata.n_obs:,}")
    print("="*80)
    
    for var in variables:
        if var not in adata.obs.columns:
            print(f"⚠️  Variable '{var}' not found in adata.obs")
            continue
            
        print(f"\n📊 {var.upper()}")
        print("-" * 50)
        
        # Get value counts
        value_counts = adata.obs[var].value_counts(dropna=False)
        percentages = adata.obs[var].value_counts(normalize=True, dropna=False) * 100
        
        # Store results
        results[var] = {
            'categories': value_counts.index.tolist(),
            'counts': value_counts.values.tolist(),
            'percentages': percentages.values.tolist(),
            'n_categories': len(value_counts),
            'missing_values': adata.obs[var].isna().sum()
        }
        
        # Display results
        summary_df = pd.DataFrame({
            'Category': value_counts.index,
            'Count': value_counts.values,
            'Percentage': [f"{p:.1f}%" for p in percentages.values]
        })
        
        print(summary_df.to_string(index=False))
        
        # Additional statistics
        print(f"\n📈 Summary:")
        print(f"   • Total categories: {len(value_counts)}")
        print(f"   • Missing values: {adata.obs[var].isna().sum()}")
        print(f"   • Most common: {value_counts.index[0]} ({value_counts.iloc[0]:,} cells, {percentages.iloc[0]:.1f}%)")
        
        if len(value_counts) > 1:
            print(f"   • Least common: {value_counts.index[-1]} ({value_counts.iloc[-1]:,} cells, {percentages.iloc[-1]:.1f}%)")
    
    return results

# Function to explore patient-level statistics
def explore_patient_level(adata, patient_id_col='Patient_ID', sample_id_col='Sample_ID'):
    """Explore patient and sample level statistics"""
    print("\n" + "="*80)
    print("PATIENT & SAMPLE LEVEL STATISTICS")
    print("="*80)
    
    # Patient level
    n_patients = adata.obs[patient_id_col].nunique()
    print(f"👥 Total patients: {n_patients}")
    
    # Sample level  
    n_samples = adata.obs[sample_id_col].nunique()
    print(f"🧪 Total samples: {n_samples}")
    
    # Cells per patient
    cells_per_patient = adata.obs.groupby(patient_id_col).size()
    print(f"\n📊 Cells per patient:")
    print(f"   • Mean: {cells_per_patient.mean():.0f}")
    print(f"   • Median: {cells_per_patient.median():.0f}")
    print(f"   • Range: {cells_per_patient.min():.0f} - {cells_per_patient.max():.0f}")
    
    # Samples per patient
    samples_per_patient = adata.obs.groupby(patient_id_col)[sample_id_col].nunique()
    print(f"\n🧪 Samples per patient:")
    print(f"   • Mean: {samples_per_patient.mean():.1f}")
    print(f"   • Median: {samples_per_patient.median():.1f}")
    print(f"   • Range: {samples_per_patient.min():.0f} - {samples_per_patient.max():.0f}")
    
    return {
        'n_patients': n_patients,
        'n_samples': n_samples,
        'cells_per_patient': cells_per_patient,
        'samples_per_patient': samples_per_patient
    }

# Set color
def color_setting():
    # Cancer Cell style color palettes
    cancer_cell_colors = {
        'major_types': {
            'T': '#E31A1C',
            'Epithelial': '#1F78B4', 
            'Plasma': '#33A02C',
            'B': '#FF7F00',
            'Myeloid': '#6A3D9A',
            'Stromal': '#B15928',
            'Encothelial': '#FB9A99',
            'Mast': '#A6CEE3'
        },
        'tissue_types': {
            'tumor': '#D62728',
            'blood': '#2CA02C',
            'normal': '#1F77B4',
            'border': '#FF7F0E',
            'Unknown': '#CCCCCC'
        },
        'treatment': {
            'Anti-PD1': '#1F77B4',
            'untreated': '#FF7F0E',
            'Anti-PD1 plus CapeOx': '#2CA02C',
            'Anti-PD1 plus Celecoxib': '#D62728',
            'Unknown': '#CCCCCC'
        },
        'microsatellite': {
            'MSI': '#E31A1C',
            'MSS': '#1F78B4',
            'Unknown': '#CCCCCC'
        },
        'response': {
            'pCR': '#2CA02C',
            'non_pCR': '#D62728',
            'Unknown': '#CCCCCC'
        },
        'treatment_stage': {
            'Pre': '#1F77B4',
            'On': '#FF7F0E', 
            'Post': '#2CA02C',
            'Unknown': '#CCCCCC'
        },
        'gender': {
            'M': '#4CAF50',
            'F': '#FF9800',
            'Unknown': '#CCCCCC'
        }
    }

    return cancer_cell_colors

# Helper function to safely handle categorical data
def save_color_maps(adata, color_dict):
    """Save color mappings to AnnData object"""
    for category, colors in color_dict.items():
        adata.uns[f'{category}_colors_map'] = colors
    print("✅ Color mappings saved to AnnData object")

# Figure 1: Overall cell type composition
def plot_cell_type_composition(adata, figsize=(8, 6), save_path=None):
    """Plot cell type composition pie chart and bar plot"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Data preparation
    cell_counts = adata.obs['Major_type'].value_counts()
    colors = [cancer_cell_colors['major_types'][ct] for ct in cell_counts.index]
    
    # Pie chart
    wedges, texts, autotexts = ax1.pie(cell_counts.values, 
                                      labels=cell_counts.index,
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90,
                                      textprops={'fontsize': 8})
    
    ax1.set_title('Cell Type Distribution\n(n=1,427,414 cells)', fontsize=8, fontweight='bold', pad=10)
    
    # Bar plot
    bars = ax2.bar(range(len(cell_counts)), cell_counts.values, color=colors)
    ax2.set_xlabel('Cell Types', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Number of Cells', fontsize=10, fontweight='bold')
    ax2.set_title('Cell Type Counts', fontsize=10, fontweight='bold')
    ax2.set_xticks(range(len(cell_counts)))
    ax2.set_xticklabels(cell_counts.index, rotation=45, ha='right', fontsize=8)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Add count labels on bars
    for bar, count in zip(bars, cell_counts.values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{count/1000:.0f}K', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    return colors

# Helper function to safely handle categorical data
def safe_fillna_categorical(series, fill_value='Unknown'):
    """Safely fill NA values in categorical series"""
    if pd.api.types.is_categorical_dtype(series):
        # Add the category if it doesn't exist
        if fill_value not in series.cat.categories:
            series = series.cat.add_categories([fill_value])
        return series.fillna(fill_value)
    else:
        return series.fillna(fill_value)

# Figure 2: Clinical characteristics overview
def plot_clinical_overview(adata, figsize=(12, 8), save_path=None):
    """Plot comprehensive clinical characteristics"""
    fig, axes = plt.subplots(2, 4, figsize=figsize)
    axes = axes.flatten()
    
    # Tissue Type
    tissue_series = safe_fillna_categorical(adata.obs['Tissue_Type'], 'Unknown')
    tissue_data = tissue_series.value_counts()
    colors_tissue = [cancer_cell_colors['tissue_types'].get(t, '#CCCCCC') for t in tissue_data.index]
    axes[0].bar(range(len(tissue_data)), tissue_data.values, color=colors_tissue)
    axes[0].set_title('Tissue Type', fontweight='bold', fontsize=8)
    axes[0].set_xticks(range(len(tissue_data)))
    axes[0].set_xticklabels(tissue_data.index, rotation=45, ha='right', fontsize=6)
    axes[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Treatment Strategy
    treat_series = safe_fillna_categorical(adata.obs['Treatment_Strategy'], 'Unknown')
    treat_data = treat_series.value_counts()
    colors_treat = [cancer_cell_colors['treatment'].get(t, '#CCCCCC') for t in treat_data.index]
    axes[1].bar(range(len(treat_data)), treat_data.values, color=colors_treat)
    axes[1].set_title('Treatment Strategy', fontweight='bold', fontsize=8)
    axes[1].set_xticks(range(len(treat_data)))
    axes[1].set_xticklabels(treat_data.index, rotation=45, ha='right', fontsize=6)
    axes[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Microsatellite Status
    msi_series = safe_fillna_categorical(adata.obs['Microsatellite_Status'], 'Unknown')
    msi_data = msi_series.value_counts()
    colors_msi = [cancer_cell_colors['microsatellite'].get(m, '#CCCCCC') for m in msi_data.index]
    axes[2].bar(range(len(msi_data)), msi_data.values, color=colors_msi)
    axes[2].set_title('Microsatellite Status', fontweight='bold', fontsize=8)
    axes[2].set_xticks(range(len(msi_data)))
    axes[2].set_xticklabels(msi_data.index, rotation=45, ha='right', fontsize=6)
    axes[2].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Treatment Response
    response_series = safe_fillna_categorical(adata.obs['Response'], 'Unknown')
    response_data = response_series.value_counts()
    colors_resp = [cancer_cell_colors['response'].get(r, '#CCCCCC') for r in response_data.index]
    axes[3].bar(range(len(response_data)), response_data.values, color=colors_resp)
    axes[3].set_title('Treatment Response', fontweight='bold', fontsize=8)
    axes[3].set_xticks(range(len(response_data)))
    axes[3].set_xticklabels(response_data.index, rotation=45, ha='right', fontsize=6)
    axes[3].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Treatment Stage
    stage_series = safe_fillna_categorical(adata.obs['Treatment_Stage'], 'Unknown')
    stage_data = stage_series.value_counts()
    colors_stage = [cancer_cell_colors['treatment_stage'].get(s, '#CCCCCC') for s in stage_data.index]
    axes[4].bar(range(len(stage_data)), stage_data.values, color=colors_stage)
    axes[4].set_title('Treatment Stage', fontweight='bold', fontsize=8)
    axes[4].set_xticks(range(len(stage_data)))
    axes[4].set_xticklabels(stage_data.index, rotation=45, ha='right', fontsize=6)
    axes[4].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Gender Distribution
    gender_series = safe_fillna_categorical(adata.obs['Gender'], 'Unknown')
    gender_data = gender_series.value_counts()
    colors_gender = [cancer_cell_colors['gender'].get(g, '#CCCCCC') for g in gender_data.index]
    axes[5].bar(range(len(gender_data)), gender_data.values, color=colors_gender)
    axes[5].set_title('Gender', fontweight='bold', fontsize=8)
    axes[5].set_xticks(range(len(gender_data)))
    axes[5].set_xticklabels(gender_data.index, fontsize=6)
    axes[5].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Tumor Stage (top categories)
    tumor_stage_series = safe_fillna_categorical(adata.obs['Tumor_stage'], 'Unknown')
    tumor_stage_data = tumor_stage_series.value_counts().head(6)
    axes[6].bar(range(len(tumor_stage_data)), tumor_stage_data.values, 
                color=plt.cm.Set3(np.linspace(0, 1, len(tumor_stage_data))))
    axes[6].set_title('Tumor Stage (Top 6)', fontweight='bold', fontsize=8)
    axes[6].set_xticks(range(len(tumor_stage_data)))
    axes[6].set_xticklabels(tumor_stage_data.index, fontsize=6)
    axes[6].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Age Distribution
    age_data = adata.obs['Age'].dropna()
    axes[7].hist(age_data, bins=20, color='#1F77B4', alpha=0.7, edgecolor='black', linewidth=0.5)
    axes[7].set_title('Age Distribution', fontweight='bold', fontsize=8)
    axes[7].set_xlabel('Age (years)', fontsize=7)
    axes[7].set_ylabel('Number of Cells', fontsize=7)
    axes[7].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# Figure 3: Cell type composition by tissue type
def plot_celltype_by_tissue(adata, figsize=(8, 5), save_path=None):
    """Plot cell type composition across tissue types"""
    # Prepare data
    tissue_series = safe_fillna_categorical(adata.obs['Tissue_Type'], 'Unknown')
    tissue_celltype = pd.crosstab(tissue_series, 
                                 adata.obs['Major_type'], normalize='index') * 100
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Stacked bar plot (percentage)
    colors = [cancer_cell_colors['major_types'][ct] for ct in tissue_celltype.columns]
    tissue_celltype.plot(kind='bar', stacked=True, ax=ax1, color=colors, 
                        legend=False, width=0.8)
    ax1.set_title('Cell Type Composition by Tissue', fontweight='bold', fontsize=8)
    ax1.set_xlabel('Tissue Type', fontsize=7, fontweight='bold')
    ax1.set_ylabel('Percentage (%)', fontsize=7, fontweight='bold')
    ax1.set_xticklabels(tissue_celltype.index, rotation=45, ha='right', fontsize=6)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
    
    # Absolute counts heatmap
    tissue_celltype_counts = pd.crosstab(tissue_series, 
                                        adata.obs['Major_type'])
    
    im = ax2.imshow(tissue_celltype_counts.values, cmap='Reds', aspect='auto')
    ax2.set_xticks(range(len(tissue_celltype_counts.columns)))
    ax2.set_yticks(range(len(tissue_celltype_counts.index)))
    ax2.set_xticklabels(tissue_celltype_counts.columns, rotation=45, ha='right', fontsize=6)
    ax2.set_yticklabels(tissue_celltype_counts.index, fontsize=6)
    ax2.set_title('Cell Counts Heatmap', fontweight='bold', fontsize=8)
    
    # Add text annotations
    for i in range(len(tissue_celltype_counts.index)):
        for j in range(len(tissue_celltype_counts.columns)):
            count = tissue_celltype_counts.iloc[i, j]
            ax2.text(j, i, f'{count/1000:.0f}K' if count > 1000 else str(count),
                    ha='center', va='center', fontsize=5,
                    color='white' if count > tissue_celltype_counts.values.max()/2 else 'black')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# Figure 4: Treatment response analysis
def plot_treatment_response_analysis(adata, figsize=(10, 6), save_path=None):
    """Plot treatment response across different variables"""
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    # Response by Treatment Strategy
    treat_series = safe_fillna_categorical(adata.obs['Treatment_Strategy'], 'Unknown')
    response_series = safe_fillna_categorical(adata.obs['Response'], 'Unknown')
    response_treat = pd.crosstab(treat_series, response_series, normalize='index') * 100
    response_treat.plot(kind='bar', ax=axes[0], 
                       color=[cancer_cell_colors['response'].get(r, '#CCCCCC') for r in response_treat.columns])
    axes[0].set_title('Response by Treatment Strategy', fontweight='bold', fontsize=8)
    axes[0].set_xticklabels(response_treat.index, rotation=45, ha='right', fontsize=6)
    axes[0].set_ylabel('Percentage (%)', fontsize=7)
    axes[0].legend(fontsize=6)
    
    # Response by Microsatellite Status
    msi_series = safe_fillna_categorical(adata.obs['Microsatellite_Status'], 'Unknown')
    response_msi = pd.crosstab(msi_series, response_series, normalize='index') * 100
    response_msi.plot(kind='bar', ax=axes[1],
                     color=[cancer_cell_colors['response'].get(r, '#CCCCCC') for r in response_msi.columns])
    axes[1].set_title('Response by MSI Status', fontweight='bold', fontsize=8)
    axes[1].set_xticklabels(response_msi.index, rotation=45, ha='right', fontsize=6)
    axes[1].set_ylabel('Percentage (%)', fontsize=7)
    axes[1].legend(fontsize=6)
    
    # Cell type composition by response
    celltype_response = pd.crosstab(response_series, adata.obs['Major_type'], normalize='index') * 100
    celltype_response.plot(kind='bar', stacked=True, ax=axes[2],
                          color=[cancer_cell_colors['major_types'][ct] for ct in celltype_response.columns])
    axes[2].set_title('Cell Types by Response', fontweight='bold', fontsize=8)
    axes[2].set_xticklabels(celltype_response.index, rotation=45, ha='right', fontsize=6)
    axes[2].set_ylabel('Percentage (%)', fontsize=7)
    axes[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=5)
    
    # Treatment stage analysis
    stage_series = safe_fillna_categorical(adata.obs['Treatment_Stage'], 'Unknown')
    stage_response = pd.crosstab(stage_series, response_series, normalize='index') * 100
    stage_response.plot(kind='bar', ax=axes[3],
                       color=[cancer_cell_colors['response'].get(r, '#CCCCCC') for r in stage_response.columns])
    axes[3].set_title('Response by Treatment Stage', fontweight='bold', fontsize=8)
    axes[3].set_xticklabels(stage_response.index, rotation=45, ha='right', fontsize=6)
    axes[3].set_ylabel('Percentage (%)', fontsize=7)
    axes[3].legend(fontsize=6)
    
    # Gender analysis
    gender_series = safe_fillna_categorical(adata.obs['Gender'], 'Unknown')
    gender_response = pd.crosstab(gender_series, response_series, normalize='index') * 100
    gender_response.plot(kind='bar', ax=axes[4],
                        color=[cancer_cell_colors['response'].get(r, '#CCCCCC') for r in gender_response.columns])
    axes[4].set_title('Response by Gender', fontweight='bold', fontsize=8)
    axes[4].set_xticklabels(gender_response.index, rotation=0, fontsize=6)
    axes[4].set_ylabel('Percentage (%)', fontsize=7)
    axes[4].legend(fontsize=6)
    
    # Age vs Response
    response_patients = adata.obs.groupby(['Patient_ID', 'Response', 'Age']).size().reset_index()
    response_patients = response_patients.dropna()
    for response in response_patients['Response'].unique():
        if response != 'Unknown' and not pd.isna(response):
            data = response_patients[response_patients['Response'] == response]['Age']
            axes[5].hist(data, alpha=0.6, label=response, bins=15,
                        color=cancer_cell_colors['response'].get(response, '#CCCCCC'))
    axes[5].set_title('Age Distribution by Response', fontweight='bold', fontsize=8)
    axes[5].set_xlabel('Age (years)', fontsize=7)
    axes[5].set_ylabel('Number of Patients', fontsize=7)
    axes[5].legend(fontsize=6)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# Figure 5: Patient and sample overview
def plot_patient_sample_overview(adata, figsize=(8, 6), save_path=None):
    """Plot patient and sample level statistics"""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # Cells per patient
    cells_per_patient = adata.obs.groupby('Patient_ID').size()
    axes[0, 0].hist(cells_per_patient, bins=20, color='#1F77B4', alpha=0.7, edgecolor='black')
    axes[0, 0].set_title('Cells per Patient Distribution', fontweight='bold', fontsize=8)
    axes[0, 0].set_xlabel('Number of Cells', fontsize=7)
    axes[0, 0].set_ylabel('Number of Patients', fontsize=7)
    axes[0, 0].axvline(cells_per_patient.median(), color='red', linestyle='--', 
                      label=f'Median: {cells_per_patient.median():.0f}')
    axes[0, 0].legend(fontsize=6)
    
    # Samples per patient
    samples_per_patient = adata.obs.groupby('Patient_ID')['Sample_ID'].nunique()
    axes[0, 1].hist(samples_per_patient, bins=range(1, samples_per_patient.max()+2), 
                   color='#2CA02C', alpha=0.7, edgecolor='black')
    axes[0, 1].set_title('Samples per Patient Distribution', fontweight='bold', fontsize=8)
    axes[0, 1].set_xlabel('Number of Samples', fontsize=7)
    axes[0, 1].set_ylabel('Number of Patients', fontsize=7)
    axes[0, 1].set_xticks(range(1, samples_per_patient.max()+1))
    
    # Study composition
    study_counts = adata.obs['study'].value_counts()
    axes[1, 0].pie(study_counts.values, labels=study_counts.index, autopct='%1.1f%%',
                  colors=plt.cm.Set3(np.linspace(0, 1, len(study_counts))))
    axes[1, 0].set_title('Study Composition', fontweight='bold', fontsize=8)
    
    # Summary statistics text
    axes[1, 1].axis('off')
    summary_text = f"""
    Dataset Summary:
    
    Total Cells: {adata.n_obs:,}
    Total Patients: {adata.obs['Patient_ID'].nunique()}
    Total Samples: {adata.obs['Sample_ID'].nunique()}
    
    Cell Types: {adata.obs['Major_type'].nunique()}
    Studies: {adata.obs['study'].nunique()}
    
    Median cells/patient: {cells_per_patient.median():.0f}
    Mean cells/patient: {cells_per_patient.mean():.0f}
    
    Treatment Strategies: {adata.obs['Treatment_Strategy'].nunique()-1}
    Response Available: {(~adata.obs['Response'].isna()).sum():,} cells
    """
    axes[1, 1].text(0.1, 0.9, summary_text, transform=axes[1, 1].transAxes,
                   fontsize=7, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()