"""
PREDICT_PLOTS UTILITY MODULE: VISUALIZE ML PREDICTION RESULTS
Created on February 28th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os, inspect, re
import numpy as np

import matplotlib.pyplot as plt
from IPython.display import Image, display

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def absolute_value(val, sizes):
    """
    This function is used in the probability proportion pie chart.

    :param val: numerical value (percentage).
    :param sizes: list of sizes for pie chart segments.
    :return: formatted string.
    """
    total = sum(sizes)
    absolute = int(val / 100. * total)
    
    return f'{absolute}\n({val:.1f}%)'


def extract_id_label(href):
    """
    Simple function to extrapolate 'id' and 'label' from href in 'prefict_df'.

    :param href: long string for 'id' and 'label' combined
    :return: 'id', and 'label' short strings
    """
    
    id_match = re.search(r'<a href="([^"]+)">', href)
    label_match = re.search(r'<a href="[^"]+">([^<]+)</a>', href)
    
    if id_match and label_match:
        return id_match.group(1), label_match.group(1)
    elif id_match:
        return id_match.group(1), None
    elif label_match:
        return None, label_match.group(1)
    else:
        return href, None


def generate_prediction_plots(predict_df, disease_dir, save_plots=True):
    """
    Generates visualization plots based on prediction results.
    
    :param predicted_df: Dataframe containing prediction results 
        ['drug', 'gene', 'predicted_interaction', 'prob'].
    :param disease_directories: dictionary with paths information.
    :param save_plots: whether to save the plots to disk (default=True).
    :return: list of paths to the generated plots.
    """
    network_directory = disease_dir['networkXGB_directory']
    disease_name_label = disease_dir['disease_name']
    date_str = disease_dir['date_string']
    plot_paths = []
    
    # generate prediction outcome proportion pie chart
    num_ones = (predict_df['predicted_interaction'] == 1).sum()
    num_zeros = (predict_df['predicted_interaction'] == 0).sum()
    plt.figure(figsize=(8, 8))
    labels = 'Predicted Interactions (1)', 'No Interactions (0)'
    sizes = [num_ones, num_zeros]
    colors = ['#ff9999', '#66b3ff']
    explode = (0.1, 0)
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct=lambda val: absolute_value(val, sizes), shadow=True, startangle=90)
    plt.axis('equal')
    plt.title('Prediction Distribution')
    
    # save pie chart if requested
    if save_plots:
        pie_chart_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_interaction_pie_chart.png')
        plt.savefig(pie_chart_path, transparent=True, dpi=300)
        plot_paths.append(pie_chart_path)
    
    plt.close()
    
    # generate prediction probability distribution histogram
    interaction_df = predict_df[predict_df['predicted_interaction'] == 1]
    if len(interaction_df) > 0:  # only generate histogram if there are positive predictions
        interaction_probs = interaction_df['prob']
        plt.figure(figsize=(10, 6))
        plt.hist(interaction_probs, bins=30, edgecolor='k', alpha=0.7, color='#4472C4')
        plt.title('Distribution of Probabilities for Predicted Interactions')
        plt.xlabel('Probability')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        
        # add mean line
        mean_prob = interaction_probs.mean()
        plt.axvline(mean_prob, color='r', linestyle='dashed', linewidth=1, 
                    label=f'Mean: {mean_prob:.3f}')
        plt.legend()
        
        # save histogram if requested
        if save_plots:
            hist_plot_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_probability_distribution.png')
            plt.savefig(hist_plot_path, transparent=True, dpi=300)
            plot_paths.append(hist_plot_path)
        
        plt.close()
    
    # create rank visualization for top drugs (if any positive predictions)
    if len(interaction_df) > 0:
        # Get top 20 drugs by interaction count
        drug_counts = interaction_df['drug'].value_counts().head(20)

        first_drug = interaction_df['drug'].iloc[0] if len(interaction_df) > 0 else ""
        
        if '<a href=' in str(first_drug):
            # Apply the extraction function to get clean drug labels
            drug_series = interaction_df['drug'].apply(lambda href: extract_id_label(href)[1] or href)
            drug_counts = drug_series.value_counts().head(20)
        else:
            drug_counts = interaction_df['drug'].value_counts().head(20)
        
        if len(drug_counts) > 0:
            plt.figure(figsize=(12, 8))
            drug_counts.plot(kind='bar', color='#70AD47')
            plt.title('Top Drugs by Number of Predicted Interactions')
            plt.xlabel('Drug')
            plt.ylabel('Number of Interactions')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            # save bar chart if requested
            if save_plots:
                bar_chart_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_top_drugs_bar_chart.png')
                plt.savefig(bar_chart_path, transparent=True, dpi=300)
                plot_paths.append(bar_chart_path)
            
            plt.close()
    
    # return the list of generated plot file paths
    return plot_paths


def print_plots(predicted_df, disease_directory):
    """
    Displays plots based on prediction results.
    
    :param predicted_df: Dataframe containing prediction results 
        ['drug', 'gene', 'predicted_interaction', 'prob'].
    :param disease_directories: dictionary with paths information.
    :return: None.
    """
    plot_paths = generate_prediction_plots(predicted_df,disease_directory)

    # iterate through paths and print each of them
    for path in plot_paths:
        display(Image(path))


def ranked_drugs_table(ranked_d_list):
    """
    Prints the table of top ranked drugs from the ML prediction.

    :param ranked_d_list: output of 'run_networkmodel()'.
    :return: None.
    """
    print(f"\nTop {len(ranked_d_list)} drugs found for repurposing:")
    for i, drug in enumerate(ranked_d_list):
        print(f"{i+1}. {drug['label']} ({drug['id']}) - {drug['count']} predicted interactions")