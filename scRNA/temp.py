import os
import pandas as pd
import gseapy as gp
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

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

def perform_survival_analysis(
    mp_scores: pd.DataFrame,
    clinical_data: pd.DataFrame,
    mp_id: str,
    time_col: str = 'OS_time',
    status_col: str = 'OS_status',
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
    median_score = score.median()
    merged_data['Group'] = ['High' if x > median_score else 'Low' for x in score]

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