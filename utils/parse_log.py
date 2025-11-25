"""
PARSE_LOGS UTILITY MODULE: parses through the iteration log
Created on March 29th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import pandas as pd
import numpy as np
import plotly.express as px
import re

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def parse_runs_log(file_path):
    """
    Parse a runs log file and return a DataFrame with extracted parameters and metrics.
    
    :param file_path: path to the log file.
    :return: DataFrame containing parsed run information.
    """
    with open(file_path, 'r') as file:
        content = file.read()
    run_sections = content.split("# # # # # # # # # # # # # # # # # # # # # # # # # # #")
    runs_data = []
    
    for section in run_sections:
        if not section.strip():
            continue  # skip empty sections
        run_data = {}
        run_type_match = re.search(r"The run using the (.*?) (with|without) negative samples", section)
        if run_type_match:
            run_data['model'] = run_type_match.group(1).strip()
            run_data['negative_samples'] = run_type_match.group(2) == 'with'
        runtime_match = re.search(r"took (.*?)\.", section)
        if runtime_match:
            run_data['runtime'] = runtime_match.group(1).strip()
        drugs_match = re.search(r"Prediction summary: (\d+) new drugs in (\d+) gene-to-drug edges", section)
        if drugs_match:
            run_data['new_drugs'] = int(drugs_match.group(1))
            run_data['gene_drug_edges'] = int(drugs_match.group(2))
        cpu_lines = re.findall(r"\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+", section)
        if cpu_lines:
            core_count = len(cpu_lines)
            run_data['cpu_cores'] = core_count
            avg_cpu_usage = sum([float(line.split()[-1]) for line in cpu_lines]) / core_count
            run_data['avg_cpu_usage'] = round(avg_cpu_usage, 2)
        total_edges_match = re.search(r"Total edges checked: (\d+)", section)
        if total_edges_match:
            run_data['total_edges_checked'] = int(total_edges_match.group(1))
        correctly_predicted_match = re.search(r"Correctly predicted: (\d+) \(([\d.]+)%\)", section)
        if correctly_predicted_match:
            run_data['correctly_predicted'] = int(correctly_predicted_match.group(1))
            run_data['success_rate'] = float(correctly_predicted_match.group(2))
        if "Parameters Used:" in section:
            param_section = section.split("Parameters Used:")[1].strip()
            param_lines = param_section.split("\n")[2:]
            for line in param_lines:
                if "-" * 5 in line or not line.strip():
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    param, value = parts
                    run_data[param] = value.strip()
        runs_data.append(run_data)
        
    df = pd.DataFrame(runs_data)
    all_params = [
        'ds_k', 'ds_min', 'ds_radius', 'ds_feat_length', 
        'ns_toggle', 'ns_method', 'ns_simil_t', 'ns_centr_t',
        'emb_jobs', 'ml_method', 'ml_model_type', 'ml_jobs',
        'ml_depth', 'ml_seed', 'ml_batch_size', 'ml_prob_filter',
        'ml_clust_filter', 'ml_iterations', 'ml_iter_max_edges'
    ]
    
    for param in all_params:
        if param not in df.columns:
            df[param] = 'na'
    
    return df


def runtime_to_seconds(runtime_str):
    """
    Convert runtime string to seconds.
    
    :param runtime_str: Runtime expressed a string of hours, minutes and seconds.
    :return: Runtime converted back to just seconds.
    """
    total_seconds = 0
    hours_match = re.search(r"(\d+) hour", runtime_str)
    if hours_match:
        total_seconds += int(hours_match.group(1)) * 3600
    minutes_match = re.search(r"(\d+) minute", runtime_str)
    if minutes_match:
        total_seconds += int(minutes_match.group(1)) * 60
    seconds_match = re.search(r"(\d+) second", runtime_str)
    if seconds_match:
        total_seconds += int(seconds_match.group(1))
    return total_seconds


def plot_runs_3d(df):
    """
    Create a 3D scatter plot of the runs with model type as color.
    
    :param df: DataFrame with parsed run information.
    :return: 3D scatterplot figure object.
    """
    plot_df = df.copy()
    plot_df['ml_iterations'] = plot_df['ml_iterations'].replace('defaulted to 1', '1')
    plot_df['ml_iterations'] = plot_df['ml_iterations'].replace('na', '1')
    plot_df['ml_iterations'] = plot_df['ml_iterations'].astype(float)
    plot_df['runtime_seconds'] = plot_df['runtime'].apply(runtime_to_seconds)
    model_colors = {'XGBoost': 'red', 'Torch-geometric (PyG)': 'blue'}
    
    fig = px.scatter_3d(
        plot_df,
        x='ml_iterations',
        y='runtime_seconds',
        z='avg_cpu_usage',
        color='model',
        color_discrete_map=model_colors,
        size='success_rate' if 'success_rate' in plot_df.columns else None,
        hover_data=['model', 'runtime', 'new_drugs', 'gene_drug_edges', 'correctly_predicted'],
        labels={
            'ml_iterations': 'Iterations',
            'runtime_seconds': 'Runtime (seconds)',
            'avg_cpu_usage': 'Avg CPU Usage (%)'
        },
        title='ML Model Performance Comparison'
    )
    
    fig.update_layout(
        scene=dict(
            xaxis_title='Iterations',
            yaxis_title='Runtime (seconds)',
            zaxis_title='Avg CPU Usage (%)',
            camera=dict(
                eye=dict(x=2, y=-2, z=1.3),
                center=dict(x=0, y=0, z=0),
                up=dict(x=0, y=0, z=1)
            )
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )
    
    return fig


def analyse_runtime(df):
    """
    Calculate average runtime and CPU usage for each model type and depth level.
    Runtime is normalized by the number of iterations.
    
    :df: DataFrame with parsed run information.
    :return: Summary of average metrics by model and depth.
    """
    analysis_df = df.copy()
    analysis_df['runtime_seconds'] = analysis_df['runtime'].apply(runtime_to_seconds)
    analysis_df['ml_depth'] = analysis_df['ml_depth'].str.replace('defaulted to ', '')
    
    for i, depth in enumerate(analysis_df['ml_depth']):
        if isinstance(depth, str) and 'ultralight' in depth:
            analysis_df.loc[i, 'ml_depth'] = 'ultralight'
        elif isinstance(depth, str) and 'light' in depth:
            analysis_df.loc[i, 'ml_depth'] = 'light'
        elif isinstance(depth, str) and 'full' in depth:
            analysis_df.loc[i, 'ml_depth'] = 'full'
            
    analysis_df['ml_iterations'] = analysis_df['ml_iterations'].replace('defaulted to 1', '1')
    analysis_df['ml_iterations'] = analysis_df['ml_iterations'].replace('na', '1')
    analysis_df['ml_iterations'] = analysis_df['ml_iterations'].astype(float)
    analysis_df['runtime_per_iteration'] = analysis_df['runtime_seconds'] / analysis_df['ml_iterations']
    summary = analysis_df.groupby(['model', 'ml_depth']).agg({
        'runtime_per_iteration': ['mean', 'min', 'max'],
        'avg_cpu_usage': ['mean', 'min', 'max']
    }).reset_index()
    summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]
    for col in ['runtime_per_iteration_mean', 'runtime_per_iteration_min', 'runtime_per_iteration_max']:
        summary[col.replace('runtime_per_iteration', 'runtime')] = summary[col].apply(
            lambda x: f"{int(x // 3600)}h {int((x % 3600) // 60)}m {int(x % 60)}s"
        )
    for col in ['avg_cpu_usage_mean', 'avg_cpu_usage_min', 'avg_cpu_usage_max']:
        summary[col] = summary[col].round(2)

    final_columns = [
        'model', 'ml_depth', 
        'runtime_mean', 'runtime_min', 'runtime_max',
        'avg_cpu_usage_mean', 'avg_cpu_usage_min', 'avg_cpu_usage_max'
    ]
    
    return summary[final_columns]
