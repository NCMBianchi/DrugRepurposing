"""
NETWORK MODEL XGB MODULE: RUNS UNSUPERVISED LEARNING WITH XGBOOST
Updated on February 27th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,multiprocessing,re,ast,warnings,datetime,inspect,gc,psutil
import numpy as np
import pandas as pd

from scipy.stats import uniform
from sklearn.model_selection import RandomizedSearchCV, RepeatedStratifiedKFold
from sklearn.utils import class_weight
from xgboost import XGBClassifier

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def format_duration(duration):
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = int(duration % 60)
    
    parts = []
    if hours > 0:
        parts.append(f"{hours} hour{'s' if hours > 1 else ''}")
    if hours > 0 or minutes > 0:
        parts.append(f"{minutes} minute{'s' if minutes > 1 else ''}")
    parts.append(f"{seconds} second{'s' if seconds > 1 else ''}")
    
    return " and ".join(parts)


def current_function_name():
    return inspect.currentframe().f_back.f_code.co_name


def unique_elements(nonUnique_list):
    """
    Short function that remove duplicate elements.
    If the list contains nodes, it will simply convert it into a set{}.
    If the list contains edges, it will remove also edges where subject and object
    are inverted, therefore not being recognised as the same by Python.

    :param nonUnique_list: biomedical entities list, where each entity is either a
        node or an edge in association networks.
    :return: list of the same biomedical entities without duplicates.
    """
    
    # if nonUnique_list is empty
    if not nonUnique_list:
        return []
    
    if isinstance(nonUnique_list[0], dict):
        # Handle list of nodes
        nodes_set = set(tuple(sorted(node.items())) for node in nonUnique_list)
        unique_list = [dict(node) for node in nodes_set]

    elif len(nonUnique_list[0]) == 4 and isinstance(nonUnique_list[0], list):
        # Handle list of edges
        unique_list = []
        seen_edges = set()
        for edge in nonUnique_list:
            subj_id = edge[0]['id']
            obj_id = edge[2]['id']
            norm_edge = tuple(sorted([subj_id, obj_id]))
            if norm_edge not in seen_edges:

                # locally store the simplified/normalised edge for parsing
                seen_edges.add(norm_edge)

                # return the actual full edge
                unique_list.append(edge)
                
    else:
        raise ValueError("Input is not recognised.")
    
    return unique_list


def extract_id_label(href):
    """
    Simple function to extrapolate 'id' and 'label' from href in 'prefict_df'.

    :param href: long string for 'id' and 'label' combined.
    :return: 'id', and 'label' short strings.
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


def check_memory():
    """
    Check available memory to determine the appropriate batch size to avoid the kernel
    dying (about 20% of the available memory).

    :return: the suggested batch size.
    """
    available_gb = psutil.virtual_memory().available / (1024**3)
    
    if available_gb > 16:
        batch_size = 5000
    elif available_gb > 8:
        batch_size = 3000
    elif available_gb > 4: 
        batch_size = 1500
    else:
        batch_size = 500
        
    print(f"Suggested batch size: {batch_size}")
    return batch_size

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def filter_rank_drugs(predict_df, prob_threshold=0.65, cluster_threshold=0.8, over_under=1):
    """
    This function filters predicted drug-gene interactions based on probability
    and ranks drugs based on their interactions with genes.
    
    :param predict_df: dataframe containing predictions with columns 
        ['drug', 'gene', 'predicted_interaction', 'prob']
    :param prob_threshold: probability threshold for filtering predictions
    :param cluster_threshold: cluster threshold for ranking drugs
    :param over_under: toggle to filter above or below cluster_threshold
    :return: ranked list of drugs with counter, and list of predicted edges
    """
    # filter the dataframe
    filt_df = predict_df[(predict_df['predicted_interaction'] == 1) & (predict_df['prob'] >= prob_threshold)]
    
    # rank drugs based on how many genes each drug interacts with
    drug_counts = filt_df['drug'].value_counts()
    
    if len(drug_counts) == 0:
        return [], []  # return empty lists if no drugs meet the criteria
        
    max_count = drug_counts.max()
    threshold_count = max_count * cluster_threshold
    
    if over_under == 1:
        ranked_drugs = drug_counts[drug_counts >= threshold_count].index.tolist()
    elif over_under == 0:
        ranked_drugs = drug_counts[drug_counts <= threshold_count].index.tolist()
    
    # create the ranked list with counts
    ranked_list = []
    
    # check if the drug values contain HTML href formatting
    first_drug = drug_counts.index[0] if len(drug_counts) > 0 else ""
    
    if '<a href=' in str(first_drug):
        # handle HTML formatted entries
        for drug_href in ranked_drugs:
            drug_id, drug_label = extract_id_label(drug_href)
            if drug_id and drug_label:
                ranked_list.append({'id': drug_id, 'label': drug_label, 'count': drug_counts[drug_href]})
    else:
        # handle plain ID entries (no HTML formatting)
        for drug_id in ranked_drugs:
            # Use the drug ID as the label if no separate label info is available
            ranked_list.append({'id': drug_id, 'label': drug_id, 'count': drug_counts[drug_id]})
    
    # create a ranked dataframe
    ranked_df = filt_df[filt_df['drug'].isin(ranked_drugs)]
    
    # build the list of new 'edges'
    predict_edges = []
    
    if '<a href=' in str(first_drug):
        # parse HTML formatted entries
        for _, row in ranked_df.iterrows():
            gene_href = row['gene']
            drug_href = row['drug']
            gene_id, gene_label = extract_id_label(gene_href)
            drug_id, drug_label = extract_id_label(drug_href)
            
            edge = [
                {'id': gene_id, 'label': gene_label},
                {'label': 'xgboost:has_association'},
                {'id': drug_id, 'label': drug_label},
                {'notes': str(row['prob'])}
            ]
            predict_edges.append(edge)
    else:
        # handle plain ID entries
        for _, row in ranked_df.iterrows():
            gene_id = row['gene']
            drug_id = row['drug']
            
            edge = [
                {'id': gene_id, 'label': gene_id},
                {'label': 'xgboost:has_association'},
                {'id': drug_id, 'label': drug_id},
                {'notes': str(row['prob'])}
            ]
            predict_edges.append(edge)
    
    return ranked_list, predict_edges


def ml_prediction(training_df, prediction_df, disease_directories, param_toggle=1, input_jobs=None,
                  depth=None, input_seed=None, batch_size=None):
    """
    This function builds a supervised learning model using the provided training
    data to predict new drug-gene interactions.
    
    :param training_df: dataframe of (gene, drug, embedding, class) used for training.
    :param prediction_df: dataframe of (gene, drug, embedding, class) used for prediction.
    :param disease_directories: base paths to where data is stored.
    :param param_toggle: toggle for loading optimized parameters.
    :param input_jobs: how many CPU cores to use.
    :param depth: parameter optimization depth ('full', 'light', 'ultralight').
    :param input_seed: random seed for reproducibility.
    :param batch_size: size of the batches for prediction.
    :return: dataframe of newly predicted interactions.
    """
    # initialize paths
    network_directory = disease_directories['networkXGB_directory']
    disease_name_label = disease_directories['disease_name']
    date_str = disease_directories['date_string']

    input_file_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_random_seed.txt')
    predict_report_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_prediction.txt')
    
    if input_seed is None or input_seed == 'random':
        input_seed = np.random.randint(1, 1e6)
        with open(input_file_path, 'w') as file:
            file.write(f"Random seed generated for the run: {input_seed}\n\n")

    if depth is None:
        depth = 'light'

    if input_jobs is None:
        input_jobs = multiprocessing.cpu_count() // 2

    print(f"The ML step is run on '{input_jobs}' cores with seed '{input_seed}'. A ({depth}) parameter optimisation step is performed, using the XGBoost method.")

    ml_train_start_time = time.time()
    
    # MODEL TRAINING
    # X = embeddings (turn list of embeddings into columns)
    emb_col = training_df['fused_embedding']
    X = pd.DataFrame(emb_col.tolist())
    X = X.astype(float)

    # create two different label transformations
    # y_1: Convert `0` to `1` and `-1` to `0` in train_df (valid negative as positive)
    y_1 = training_df['class'].replace({0: 1, -1: 0}).astype(int)
    
    # y_2: Convert only `-1` to `0` in train_df (valid negative as negative)
    y_2 = training_df['class'].replace({-1: 0}).astype(int)
    
    # define parameters based on depth
    if depth == 'full':
        parameters = {
            'min_child_weight': [2, 3, 5, 8, 13, 20, 30],
            'gamma': [0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.5, 4, 6],
            'reg_alpha': [0, 0.5, 1, 3, 5, 10],
            'reg_lambda': [0, 0.5, 1, 3, 5, 10],
            'subsample': uniform(0.5, 0.5),
            'colsample_bytree': uniform(0.2, 0.8),
            'max_depth': [4, 6, 8, 10, 12, 14, 16],
            'n_estimators': [35, 45, 50, 70, 80, 90, 100],
            'learning_rate': uniform(0, 0.3),
        }
        n_s = 10
        n_r = 5
        n_i = 20
    elif depth == 'light':
        parameters = {
            'min_child_weight': [3, 8, 20],
            'gamma': [0.5, 1.2, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20, 30],
            'learning_rate': [0.1],
        }
        n_s = 3
        n_r = 2
        n_i = 5
    elif depth == 'ultralight':
        parameters = {
            'min_child_weight': [5, 20],
            'gamma': [0.5, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20],
            'learning_rate': [0.1],
        }
        n_s = 2
        n_r = 1
        n_i = 5

    params_file_path_1 = os.path.join(network_directory, f'{disease_name_label}_{date_str}_hyperparameters_1.txt')
    params_file_path_2 = os.path.join(network_directory, f'{disease_name_label}_{date_str}_hyperparameters_2.txt')

    # check for existing parameters file
    if param_toggle == 1 and os.path.exists(params_file_path_1):
        try:
            with open(params_file_path_1, 'r') as file:
                file_content = file.read().strip()
                if file_content:  # Check if file is not empty
                    param_dict = ast.literal_eval(file_content)
                    parameters = {k: [v] for k, v in param_dict.items()}
                    n_s = 2
                    n_r = 1
                    n_i = 1
                    print("Loaded parameters from existing file.")
                else:
                    print("Parameter file exists but is empty, using default parameters.")
        except Exception as e:
            print(f"Error loading parameters file: {str(e)}")
            print("Using default parameters instead.")

    # train model_1: Valid negative (0) -> positive (1), Invalid negative (-1) -> negative (0)
    xgb_model_1 = XGBClassifier(objective='multi:softmax', eval_metric='mlogloss', 
                                num_class=2, random_state=input_seed, n_jobs=input_jobs,
                               tree_method='hist')
    rskf_1 = RepeatedStratifiedKFold(n_splits=n_s, n_repeats=n_r, random_state=input_seed)
    randomized_search_1 = RandomizedSearchCV(xgb_model_1, param_distributions=parameters,
                                            scoring='f1_weighted', n_iter=n_i, n_jobs=input_jobs,
                                            error_score='raise', cv=rskf_1.split(X, y_1), verbose=1,
                                            refit=True, random_state=input_seed)
    weight_1 = class_weight.compute_sample_weight('balanced', y_1)
    randomized_search_1.fit(X, y_1, sample_weight=weight_1)
    best_model_1 = randomized_search_1.best_estimator_

    del rskf_1, xgb_model_1
    gc.collect()

    # train model_2: Both valid (0) and invalid (-1) negative -> negative (0)
    xgb_model_2 = XGBClassifier(objective='multi:softmax', eval_metric='mlogloss', 
                               num_class=2, random_state=input_seed, n_jobs=input_jobs,
                               tree_method='hist')
    rskf_2 = RepeatedStratifiedKFold(n_splits=n_s, n_repeats=n_r, random_state=input_seed)
    randomized_search_2 = RandomizedSearchCV(xgb_model_2, param_distributions=parameters,
                                           scoring='f1_weighted', n_iter=n_i, n_jobs=input_jobs,
                                           error_score='raise', cv=rskf_2.split(X, y_2), verbose=1,
                                           refit=True, random_state=input_seed)
    weight_2 = class_weight.compute_sample_weight('balanced', y_2)
    randomized_search_2.fit(X, y_2, sample_weight=weight_2)
    best_model_2 = randomized_search_2.best_estimator_

    del rskf_2, xgb_model_2
    gc.collect()

    # save best parameters
    best_params_1 = randomized_search_1.best_params_
    best_score_1 = randomized_search_1.best_score_
    best_params_2 = randomized_search_2.best_params_
    best_score_2 = randomized_search_2.best_score_
    with open(params_file_path_1, 'w') as file:
        file.write(str(best_params_1))
    with open(params_file_path_2, 'w') as file:
        file.write(str(best_params_2))

    del randomized_search_1, randomized_search_2
    del X, y_1, y_2, parameters, weight_1, weight_2
    gc.collect()

    ml_train_end_time = time.time()
    train_duration = ml_train_end_time - ml_train_start_time
    train_formatted_duration = format_duration(train_duration)

    # PREDICTION
    ml_predict_start_time = time.time()
    with open(predict_report_path, 'w') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting prediction process at {current_time}...\n\n")

    # process in smaller batches
    if batch_size == None:
        batch_size = check_memory()
    else:
        print(f"Input batch size: {batch_size}")

    total_rows = len(prediction_df)
    total_batches = (total_rows + batch_size - 1) // batch_size
    final_predictions = []
    final_probabilities = []

    report_interval = max(1, total_batches // 10)
    last_report = -1

    for i in range(0, total_rows, batch_size):
        batch_end = min(i + batch_size, total_rows)
        batch_num = i // batch_size
        
        # extract and convert just this batch of embeddings
        batch_df = prediction_df.iloc[i:batch_end]
        emb_batch = batch_df['fused_embedding']
        X_batch = pd.DataFrame(emb_batch.tolist()).astype(float)
        
        # predict with both models
        predictions_1 = best_model_1.predict(X_batch)
        predictions_prob_1 = best_model_1.predict_proba(X_batch)
        
        predictions_2 = best_model_2.predict(X_batch)
        predictions_prob_2 = best_model_2.predict_proba(X_batch)
        
        # combine predictions
        batch_predictions = np.logical_or(predictions_1, predictions_2).astype(int)
        batch_probabilities = np.max((predictions_prob_1 + predictions_prob_2) / 2, axis=1)
        
        # store results
        final_predictions.extend(batch_predictions)
        final_probabilities.extend(batch_probabilities)
        
        # force garbage collection
        del X_batch, predictions_1, predictions_prob_1, predictions_2, predictions_prob_2
        gc.collect()
        
        # log progress at 10% intervals
        current_percentage = batch_num * 100 // total_batches
        if current_percentage // 10 > last_report:
            last_report = current_percentage // 10
            with open(predict_report_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"Progress: {current_percentage}% ({batch_num+1}/{total_batches} batches) processed at {current_time}\n")
    
    # create the final dataframe
    interaction_predictions_df = pd.DataFrame({
        'drug': prediction_df['drug'],
        'gene': prediction_df['gene'],
        'predicted_interaction': final_predictions,
        'prob': final_probabilities
    })

    # save predictions
    interaction_predictions_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_predicted_df.csv')
    interaction_predictions_df.to_csv(interaction_predictions_path, index=False)

    ml_predict_end_time = time.time()
    predict_duration = ml_predict_end_time - ml_predict_start_time
    predict_formatted_duration = format_duration(predict_duration)

    run_times = [train_formatted_duration, predict_formatted_duration]

    return interaction_predictions_df, run_times
    

def run_network_model_xgb(gene_embeddings, drug_embeddings, alldrug_embeddings, 
                     training_df, prediction_df, nodes, drug_nodes, disease_directories,
                     ns_toggle=1, run_jobs=None, run_depth=None,
                     run_seed=None, prob_input=None, clust_input=None, batch_s=None,
                     param_toggle=1, ou_toggle=1):
    """
    This function runs the entire network_model script and saves nodes and edges files
    after a prediction made via the chosen ML method.

    :param gene_embeddings: gene embeddings from run_embeddings().
    :param drug_embeddings: disease drug embeddings from run_embeddings().
    :param alldrug_embeddings: all drug embeddings from run_embeddings().
    :param training_df: training dataframe from run_embeddings().
    :param prediction_df: prediction dataframe from run_embeddings().
    :param nodes: all nodes from previous steps.
    :param drug_nodes: all drug nodes from DGIdb.
    :param disease_directories: base paths to where data is stored.
    :param ns_toggle: toggle for negative samples (1 if used, 0 if not).
    :param run_jobs: how many CPU cores to use.
    :param run_depth: parameter optimization depth.
    :param run_seed: random seed for reproducibility.
    :param prob_input: probability threshold for filtering predictions.
    :param clust_input: cluster threshold for ranking drugs.
    :param batch_s: size of the batches for prediction.
    :param param_toggle: toggle for loading parameter files.
    :param ou_toggle: toggle for filtering above or below cluster threshold.
    :return: prediction edges, nodes, ranked drugs, prediction DataFrame for subsequent plots.
    """
    start_time = time.time()
    
    print(f"NOW RUNNING: {current_function_name()} after 'run_embeddingsXGB()'.")
    
    # initialize paths
    network_directory = disease_directories['networkXGB_directory']
    disease_name_label = disease_directories['disease_name']
    date_str = disease_directories['date_string']

    # train ML model and make predictions
    predicted_df, ml_run_times = ml_prediction(training_df, prediction_df, disease_directories,
                                               param_toggle=param_toggle, input_jobs=run_jobs,
                                               depth=run_depth, input_seed=run_seed,
                                               batch_size=batch_s)
    
    training_formatted_duration = ml_run_times[0]
    prediction_formatted_duration = ml_run_times[1]

    # rank and filter predicted drug interactions
    rank_start_time = time.time()

    all_nodes = nodes + drug_nodes
    id_to_label = {node['id']: node['label'] for node in all_nodes}

    predicted_df['drug'] = predicted_df['drug'].apply(
        lambda drug_uri: f'<a href="{drug_uri}">{id_to_label.get(drug_uri, drug_uri)}</a>'
    )
    predicted_df['gene'] = predicted_df['gene'].apply(
        lambda gene_uri: f'<a href="{gene_uri}">{id_to_label.get(gene_uri, gene_uri)}</a>'
    )
    
    if not prob_input:
        prob_input = 0.65
    if not clust_input:
        clust_input = 0.8
        
    ranked_drug_list, predicted_edges = filter_rank_drugs(
        predicted_df, 
        prob_threshold=prob_input,
        cluster_threshold=clust_input, 
        over_under=ou_toggle
    )
    
    for ranked_drug_elem in ranked_drug_list:
        print(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")

    rank_end_time = time.time()
    rank_duration = rank_end_time - rank_start_time
    formatted_rank_duration = format_duration(rank_duration)
    
    # save predictions
    full_edges = predicted_edges
    full_nodes = [edge[0] for edge in full_edges] + [edge[2] for edge in full_edges]
    full_nodes = unique_elements(full_nodes)

    full_edges_df = pd.DataFrame([
        {
            'subject_id': edge[0]['id'],
            'subject_label': edge[0]['label'],
            'relation': edge[1]['label'],
            'object_id': edge[2]['id'],
            'object_label': edge[2]['label'],
            'notes': edge[3]['notes']
        }
        for edge in full_edges
    ])
    full_edges_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_prediction_edges.csv'), index=False)
    
    full_nodes_df = pd.DataFrame(full_nodes)
    full_nodes_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_prediction_nodes.csv'), index=False)

    end_time = time.time()
    duration = end_time - start_time
    formatted_duration = format_duration(duration)
    print(f"'networkmodel.py' run finished in {formatted_duration} –where the training itself took {training_formatted_duration}, the prediction took {prediction_formatted_duration}, and the ranking of drugs took {formatted_rank_duration}.")
    
    return full_edges, full_nodes, ranked_drug_list, predicted_df