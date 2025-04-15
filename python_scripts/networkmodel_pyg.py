"""
NETWORK MODEL PYG MODULE: RUNS GRAPH NEURAL NETWORKS WITH PYTORCH GEOMETRIC
Created on March 3rd 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,multiprocessing,re,ast,warnings,datetime,inspect,gc,psutil
import pandas as pd
import numpy as np

import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, SAGEConv
from torch_geometric.data import Data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

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
    if not nonUnique_list:
        return []
    
    if isinstance(nonUnique_list[0], dict):
        # handles list of nodes
        nodes_set = set(tuple(sorted(node.items())) for node in nonUnique_list)
        unique_list = [dict(node) for node in nodes_set]

    elif len(nonUnique_list[0]) == 4 and isinstance(nonUnique_list[0], list):
        # handles list of edges
        unique_list = []
        seen_edges = set()
        for edge in nonUnique_list:
            subj_id = edge[0]['id']
            obj_id = edge[2]['id']
            norm_edge = tuple(sorted([subj_id, obj_id]))
            if norm_edge not in seen_edges:
                seen_edges.add(norm_edge)
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
    import re
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class GCNLinkPredictor(torch.nn.Module):
    """
    Graph Convolutional Network for link prediction.

    :return: GCN link predictor.
    """
    def __init__(self, in_channels, hidden_channels, out_channels, dropout=0.5):
        super(GCNLinkPredictor, self).__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)
        self.link_predictor = torch.nn.Sequential(
            torch.nn.Linear(out_channels * 2, hidden_channels),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden_channels, 1)
        )
        
    def encode(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)
        return x
        
    def decode(self, z, edge_index):
        row, col = edge_index
        return self.predict_link(z[row], z[col])
        
    def predict_link(self, z_src, z_dst):
        # Combine node embeddings for link prediction
        z_combined = torch.cat([z_src, z_dst], dim=1)
        return torch.sigmoid(self.link_predictor(z_combined))
    
    def forward(self, x, edge_index, pred_edge_index):
        # Get node embeddings
        z = self.encode(x, edge_index)
        
        # Predict links
        link_pred = self.decode(z, pred_edge_index)
        
        return link_pred


class GATLinkPredictor(torch.nn.Module):
    """
    Graph Attention Network for link prediction.

    :return: GAT link predictor.
    """
    def __init__(self, in_channels, hidden_channels, out_channels, heads=2, dropout=0.5):
        super(GATLinkPredictor, self).__init__()
        self.conv1 = GATConv(in_channels, hidden_channels, heads=heads, dropout=dropout)
        self.conv2 = GATConv(hidden_channels * heads, out_channels, heads=1, concat=False, dropout=dropout)
        self.link_predictor = torch.nn.Sequential(
            torch.nn.Linear(out_channels * 2, hidden_channels),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden_channels, 1)
        )
        
    def encode(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.elu(x)
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)
        return x
        
    def decode(self, z, edge_index):
        row, col = edge_index
        return self.predict_link(z[row], z[col])
        
    def predict_link(self, z_src, z_dst):
        z_combined = torch.cat([z_src, z_dst], dim=1)
        return torch.sigmoid(self.link_predictor(z_combined))
    
    def forward(self, x, edge_index, pred_edge_index):
        z = self.encode(x, edge_index)
        link_pred = self.decode(z, pred_edge_index)
        return link_pred


class SAGELinkPredictor(torch.nn.Module):
    """
    GraphSAGE for link prediction.

    :return: SAGE link predictor.
    """
    def __init__(self, in_channels, hidden_channels, out_channels, dropout=0.5):
        super(SAGELinkPredictor, self).__init__()
        self.conv1 = SAGEConv(in_channels, hidden_channels)
        self.conv2 = SAGEConv(hidden_channels, out_channels)
        self.link_predictor = torch.nn.Sequential(
            torch.nn.Linear(out_channels * 2, hidden_channels),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(hidden_channels, 1)
        )
        
    def encode(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv2(x, edge_index)
        return x
        
    def decode(self, z, edge_index):
        row, col = edge_index
        return self.predict_link(z[row], z[col])
        
    def predict_link(self, z_src, z_dst):
        z_combined = torch.cat([z_src, z_dst], dim=1)
        return torch.sigmoid(self.link_predictor(z_combined))
    
    def forward(self, x, edge_index, pred_edge_index):
        z = self.encode(x, edge_index)
        link_pred = self.decode(z, pred_edge_index)
        return link_pred


def train_model(model, train_data, optimizer, device):
    """
    Train the PyG model.
    
    :param model: PyG model to train.
    :param train_data: training data with graph structure and edge labels.
    :param optimizer: PyTorch optimizer.
    :param device: device to run the model on (CPU or GPU).
    :return: training loss.
    """
    model.train()
    optimizer.zero_grad()
    
    # move data to device
    x = train_data.x.to(device)
    edge_index = train_data.edge_index.to(device)
    train_edge_index = train_data.train_edge_index.to(device)
    train_edge_label = train_data.train_edge_label.to(device)
    
    # forward pass
    pred = model(x, edge_index, train_edge_index).squeeze()
    
    # calculate loss
    loss = F.binary_cross_entropy(pred, train_edge_label)
    
    # backward pass
    loss.backward()
    optimizer.step()
    
    return loss.item()


def validate_model(model, train_data, val_edge_index, val_edge_label, device):
    """
    Validate the PyG model.
    
    :param model: PyG model to validate
    :param train_data: Training data with graph structure
    :param val_edge_index: Validation edge indices
    :param val_edge_label: Validation edge labels
    :param device: Device to run the model on (CPU or GPU)
    :return: Validation loss
    """
    model.eval()
    
    # move data to device
    x = train_data.x.to(device)
    edge_index = train_data.edge_index.to(device)
    val_edge_index = val_edge_index.to(device)
    val_edge_label = val_edge_label.to(device)
    
    with torch.no_grad():
        # forward pass
        pred = model(x, edge_index, val_edge_index).squeeze()
        
        # calculate loss
        loss = F.binary_cross_entropy(pred, val_edge_label)
    
    return loss.item()


def predict_links(model, train_data, pred_data, device, batch_size=5000, r_path=None):
    """
    Predicts links using the trained model, filtering out any pairs with indices
    that are out of bounds for the training data.
    
    :param model: trained PyG model (GCN, GAT, or SAGE).
    :param train_data: training data with graph structure and node features.
    :param pred_data: prediction data with edge indices to predict.
    :param device: device to run the model on (CPU or GPU).
    :param batch_size: batch size for prediction to avoid memory issues.
    :param r_path: path to write progress reports.
    :return: tuple of prediction probabilities, and prediction binary label.
    """
    model.eval()
    
    # get prediction edge index
    pred_edge_index = pred_data.pred_edge_index
    
    # get total number of edges to predict
    total_edges = pred_edge_index.size(1)
    total_batches = (total_edges + batch_size - 1) // batch_size
    
    # calculate reporting interval
    report_interval = max(1, total_batches // 10)
    
    # initialize predictions
    all_preds = []
    
    # get size of z tensor to check bounds
    z_size = train_data.x.size(0)
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Starting prediction on {total_edges} pairs in {total_batches} batches at {current_time}...\n")
    
    # filter the edge index to only include valid indices
    valid_pairs = []
    valid_indices = []
    
    # first scan to find all valid pairs
    for i in range(total_edges):
        row, col = pred_edge_index[0, i].item(), pred_edge_index[1, i].item()
        if row < z_size and col < z_size:
            valid_pairs.append((row, col))
            valid_indices.append(i)
    
    valid_count = len(valid_pairs)
    filtered_batches = (valid_count + batch_size - 1) // batch_size
    
    print(f"Filtering {total_edges} pairs to {valid_count} valid pairs within bounds.")
    
    # predict in batches to avoid OOM errors
    for batch_idx in range(filtered_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, valid_count)
        
        batch_indices = valid_indices[start_idx:end_idx]
        batch_rows = [pred_edge_index[0, i].item() for i in batch_indices]
        batch_cols = [pred_edge_index[1, i].item() for i in batch_indices]
        
        batch_edge_index = torch.tensor([batch_rows, batch_cols], dtype=torch.long)
        
        # move data to device
        x = train_data.x.to(device)
        edge_index = train_data.edge_index.to(device)
        batch_edge_index = batch_edge_index.to(device)
        
        with torch.no_grad():
            # forward pass
            batch_pred = model(x, edge_index, batch_edge_index).squeeze()
            all_preds.append(batch_pred.cpu())
        
        # report progress every 10%
        if r_path and (batch_idx + 1) % report_interval == 0:
            progress_percent = min(100, (batch_idx + 1) * 100 // filtered_batches)
            with open(r_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"Prediction progress: {progress_percent}% ({batch_idx+1}/{filtered_batches} batches) processed at {current_time}\n")
    
    # concatenate all predictions
    if all_preds:
        all_preds = torch.cat(all_preds)
        
        # create full-sized prediction arrays, with zeros for filtered-out pairs
        full_probabilities = np.zeros(total_edges)
        full_predictions = np.zeros(total_edges)
        
        # fill in the predictions for valid pairs
        for i, idx in enumerate(valid_indices):
            if i < len(all_preds):
                prob = all_preds[i].item()
                full_probabilities[idx] = prob
                full_predictions[idx] = 1 if prob > 0.5 else 0
    else:
        full_probabilities = np.zeros(total_edges)
        full_predictions = np.zeros(total_edges)
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Prediction completed at {current_time}\n\n")
    
    return full_probabilities, full_predictions


def filter_rank_drugs(gene_drug_pairs, predictions, probabilities, node_labels, 
                      prob_threshold=0.65, cluster_threshold=0.8, over_under=1, r_path=None):
    """
    Filter and rank drugs based on prediction results.
    
    :param gene_drug_pairs: list of gene-drug pairs.
    :param predictions: binary prediction results.
    :param probabilities: probability prediction results.
    :param node_labels: dictionary mapping node IDs to their labels.
    :param prob_threshold: probability threshold for filtering.
    :param cluster_threshold: threshold for clustering drugs.
    :param over_under: toggle for filtering mode.
    :param r_path: path to write progress reports.
    :return: ranked drug list and predicted edges.
    """
    start_time = time.time()
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Starting filter_rank_drugs with {len(gene_drug_pairs)} pairs at {current_time}...\n")
    
    # create DataFrame for predictions
    results = []
    total_pairs = min(len(gene_drug_pairs), len(predictions))
    report_interval = max(1, total_pairs // 10)
    
    for i, pair in enumerate(gene_drug_pairs):
        if i < len(predictions):
            results.append({
                'gene_id': pair['gene_id'],
                'drug_id': pair['drug_id'],
                'predicted_interaction': predictions[i],
                'prob': probabilities[i]
            })
            
            # Report progress every 10%
            if r_path and (i + 1) % report_interval == 0:
                progress_percent = (i + 1) * 100 // total_pairs
                with open(r_path, 'a') as f:
                    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    elapsed = format_duration(time.time() - start_time)
                    f.write(f"DataFrame creation: {progress_percent}% ({i+1}/{total_pairs}), elapsed time: {elapsed} at {current_time}\n")
    
    pred_df = pd.DataFrame(results)

    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"DataFrame creation with {len(pred_df)} rows ended at {current_time}...\n")
    
    # filter predictions above threshold
    filtered_df = pred_df[(pred_df['predicted_interaction'] == 1) & (pred_df['prob'] >= prob_threshold)]
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Filtered to {len(filtered_df)} positive predictions with probability >= {prob_threshold} at {current_time}...\n")
    
    # count interactions per drug
    drug_counts = filtered_df['drug_id'].value_counts()
    
    if len(drug_counts) == 0:
        if r_path:
            with open(r_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"No drugs meet criteria, returning empty results at {current_time}\n")
        end_time = time.time()
        run_time = end_time - start_time
        formatted_run_time = format_duration(run_time)
        return [], [], formatted_run_time
    
    # get max count for threshold calculation
    max_count = drug_counts.max()
    threshold_count = max_count * cluster_threshold
    
    # filter drugs based on interaction count
    if over_under == 1:
        ranked_drugs = drug_counts[drug_counts >= threshold_count].index.tolist()
    else:
        ranked_drugs = drug_counts[drug_counts <= threshold_count].index.tolist()

    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Selected {len(ranked_drugs)} ranked drugs at {current_time}...\n")
    
    # create ranked drug list
    ranked_list = []
    for drug_id in ranked_drugs:
        # get drug label from our lookup dictionary
        drug_label = node_labels.get(drug_id, drug_id)
        
        ranked_list.append({
            'id': drug_id,
            'label': drug_label,
            'count': drug_counts[drug_id]
        })
    
    # create edges for predicted interactions
    filtered_subset = filtered_df[filtered_df['drug_id'].isin(ranked_drugs)]
    total_edges = len(filtered_subset)
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Creating {total_edges} predicted edges at {current_time}...\n")
    
    edge_report_interval = max(1, total_edges // 10)
    predicted_edges = []
    
    for i, (_, row) in enumerate(filtered_subset.iterrows()):
        gene_id = row['gene_id']
        drug_id = row['drug_id']
        
        # get labels directly from lookup dictionary
        gene_label = node_labels.get(gene_id, gene_id)
        drug_label = node_labels.get(drug_id, drug_id)
        
        edge = [
            {'id': gene_id, 'label': gene_label},
            {'label': 'pyg:has_association'},
            {'id': drug_id, 'label': drug_label},
            {'notes': str(row['prob'])}
        ]
        
        predicted_edges.append(edge)
        
        # report progress every 10%
        if r_path and (i + 1) % edge_report_interval == 0:
            progress_percent = (i + 1) * 100 // total_edges
            with open(r_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                elapsed = format_duration(time.time() - start_time)
                f.write(f"Edge creation: {progress_percent}% ({i+1}/{total_edges}), total elapsed: {elapsed} at {current_time}\n")
    
    if r_path:
        with open(r_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Entire filter_rank_drugs() function with {len(predicted_edges)} predicted edges ended at {current_time}\n\n")

    end_time = time.time()
    run_time = end_time - start_time
    formatted_run_time = format_duration(run_time)
    
    return ranked_list, predicted_edges, formatted_run_time


def run_network_model_pyg(gene_embeddings, drug_embeddings, alldrug_embeddings, 
                      training_df, prediction_df, nodes, drug_nodes, disease_directories,
                      training_data, prediction_data, 
                      ns_toggle=1, model_type="gcn", run_jobs=None, run_depth=None,
                      run_seed=None, prob_input=None, clust_input=None, batch_s=2000,
                      param_toggle=1, ou_toggle=1):
    """
    Run the PyG graph neural network model for drug repurposing prediction.
    
    :param gene_embeddings: gene embeddings from run_embeddings_pyg().
    :param drug_embeddings: drug embeddings from run_embeddings_pyg().
    :param alldrug_embeddings: all drug embeddings from run_embeddings_pyg().
    :param training_df: training DataFrame (unused in PyG, kept for interface compatibility).
    :param prediction_df: prediction DataFrame (unused in PyG, kept for interface compatibility).
    :param nodes: all nodes from previous steps.
    :param drug_nodes: all drug nodes from DGIdb.
    :param disease_directories: base paths to where data is stored.
    :param training_data: PyG training data.
    :param prediction_data: PyG prediction data.
    :param ns_toggle: toggle for negative samples.
    :param model_type: type of GNN model to use.
    :param run_jobs: CPU cores for training.
    :param run_depth: training depth.
    :param run_seed: random seed.
    :param prob_input: probability threshold for filtering.
    :param clust_input: cluster threshold for ranking.
    :param batch_s: batch size for prediction.
    :param param_toggle: toggle for loading parameters.
    :param ou_toggle: toggle for threshold mode.
    :return: prediction edges, nodes, ranked drugs, and prediction DataFrame.
    """
    start_time = time.time()
    
    print(f"NOW RUNNING: {current_function_name()} after 'run_embeddingsPYG()'.")
    
    # initialize paths
    network_directory = disease_directories['networkPYG_directory']
    disease_name_label = disease_directories['disease_name']
    date_str = disease_directories['date_string']
    train_report_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_training.txt')
    pred_report_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_prediction.txt')
    model_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_{model_type}_model.pt')

    # node label dictionary setup
    node_labels = {}
    for node_id, label in prediction_data.node_labels.items():
        node_labels[node_id] = label
    print(f"Created node label dictionary with {len(node_labels)} entries.")
    
    # set probability and cluster thresholds
    if not prob_input:
        prob_input = 0.65
    if not clust_input:
        clust_input = 0.8
    
    # set random seed
    if run_seed:
        torch.manual_seed(run_seed)
        np.random.seed(run_seed)
    
    # choose device (CPU or GPU)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    
    # check if we should load an existing model
    should_train = True
    if param_toggle == 1 and os.path.exists(model_path):
        should_train = False
        print(f"Loading existing model from {model_path}")
    
    # define training parameters based on depth
    if run_depth == "full":
        epochs = 200
        hidden_channels = 128
        patience = 20
    elif run_depth == "light":
        epochs = 100
        hidden_channels = 64
        patience = 10
    else:  # ultralight
        epochs = 50
        hidden_channels = 32
        patience = 5
    
    # get input and output dimensions
    in_channels = training_data.x.size(1)
    out_channels = 64
    
    # initialize model
    if model_type == "gcn":
        model = GCNLinkPredictor(in_channels, hidden_channels, out_channels)
    elif model_type == "gat":
        model = GATLinkPredictor(in_channels, hidden_channels, out_channels)
    elif model_type == "sage":
        model = SAGELinkPredictor(in_channels, hidden_channels, out_channels)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
    
    # move model to device
    model = model.to(device)
    
    # initialize optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    
    # training phase
    if should_train:
        train_start_time = time.time()

        with open(train_report_path, 'w') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Starting model training at {current_time}...\n")
        
        # split training edges for validation
        num_train_edges = training_data.train_edge_index.size(1)
        indices = torch.randperm(num_train_edges)
        
        # use 80% for training, 20% for validation
        train_size = int(0.8 * num_train_edges)
        train_indices = indices[:train_size]
        val_indices = indices[train_size:]
        
        # create training and validation edge indices and labels
        train_edge_index = training_data.train_edge_index[:, train_indices]
        train_edge_label = training_data.train_edge_label[train_indices]
        
        val_edge_index = training_data.train_edge_index[:, val_indices]
        val_edge_label = training_data.train_edge_label[val_indices]
        
        # store original train edge index and label
        original_train_edge_index = training_data.train_edge_index
        original_train_edge_label = training_data.train_edge_label
        
        # replace with split training data
        training_data.train_edge_index = train_edge_index
        training_data.train_edge_label = train_edge_label
        
        # initialize training log
        with open(train_report_path, 'a') as f:
            f.write(f"Training PyG {model_type.upper()} model\n")
            f.write(f"Date: {datetime.datetime.now()}\n")
            f.write(f"Epochs: {epochs}\n")
            f.write(f"Hidden channels: {hidden_channels}\n")
            f.write(f"Device: {device}\n\n")
        
        # training loop
        best_val_loss = float('inf')
        best_epoch = 0
        patience_counter = 0
        
        for epoch in range(epochs):
            # train
            loss = train_model(model, training_data, optimizer, device)
            
            # validate
            val_loss = validate_model(model, training_data, val_edge_index, val_edge_label, device)
            
            # log progress
            log_message = f"Epoch {epoch+1}/{epochs}, Train Loss: {loss:.4f}, Val Loss: {val_loss:.4f}"
            print(log_message)
            
            with open(train_report_path, 'a') as f:
                f.write(log_message + "\n")
            
            # check for improvement
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                best_epoch = epoch
                patience_counter = 0
                
                # save best model
                torch.save(model.state_dict(), model_path)
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break
        
        # restore original train edge index and label
        training_data.train_edge_index = original_train_edge_index
        training_data.train_edge_label = original_train_edge_label
        
        train_end_time = time.time()
        train_duration = train_end_time - train_start_time
        formatted_train_duration = format_duration(train_duration)
        print(f"Training completed, with best epoch: {best_epoch+1}")

        with open(train_report_path, 'a') as f:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Training completed in {formatted_train_duration} at {current_time}\n\n")
    else:
        formatted_train_duration = "0 seconds"
    
    # load best model
    model.load_state_dict(torch.load(model_path))
    
    # prediction phase
    predict_start_time = time.time()

    with open(pred_report_path, 'w') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting prediction phase at {current_time}...\n")
    
    # get gene-drug pairs for prediction
    gene_drug_pairs = prediction_data.gene_drug_pairs
    
    # predict links
    print(f"Predicting links for {len(gene_drug_pairs)} gene-drug pairs...")
    probabilities, predictions = predict_links(model, training_data, prediction_data, device, batch_s, r_path=pred_report_path)
    
   # create output for compatibility with the rest of the pipeline
    ranked_drug_list, predicted_edges, formatted_filtering_duration = filter_rank_drugs(
        gene_drug_pairs, predictions, probabilities, node_labels,
        prob_threshold=prob_input, cluster_threshold=clust_input, over_under=ou_toggle, r_path=pred_report_path
    )
    
    # save all ranked drugs to file and print limited number to console
    with open(os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_ranked_drugs.txt'), 'w') as f:
        f.write(f"Found {len(ranked_drug_list)} drugs ranked by predicted interactions:\n\n")
        for ranked_drug_elem in ranked_drug_list:
            f.write(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']} ({ranked_drug_elem['label']})\n")
    
    # print a maximum of 10 drugs with their counter
    if len(ranked_drug_list) <= 10:
        # if 10 or fewer, print all
        for ranked_drug_elem in ranked_drug_list:
            print(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")
    else:
        # otherwise print top 5 and bottom 5
        print("Top 5 ranked drugs:")
        for ranked_drug_elem in ranked_drug_list[:5]:
            print(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")
        print("(...)")
        print("Bottom 5 ranked drugs:")
        for ranked_drug_elem in ranked_drug_list[-5:]:
            print(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")
    
    # extract nodes from edges
    predicted_nodes = [edge[0] for edge in predicted_edges] + [edge[2] for edge in predicted_edges]
    predicted_nodes = unique_elements(predicted_nodes)
    
    # create prediction DataFrame for compatibility
    results = []
    for i, pair in enumerate(gene_drug_pairs):
        if i < len(predictions):
            # format for compatibility with visualization
            gene_id, drug_id = pair['gene_id'], pair['drug_id']
            
            # get labels from node_labels dictionary
            gene_label = node_labels.get(gene_id, gene_id)
            drug_label = node_labels.get(drug_id, drug_id)
            
            # create HTML format for consistency
            gene_href = f'<a href="{gene_id}">{gene_label}</a>'
            drug_href = f'<a href="{drug_id}">{drug_label}</a>'
            
            results.append({
                'gene': gene_href,
                'drug': drug_href,
                'predicted_interaction': int(predictions[i]),
                'prob': float(probabilities[i])
            })
    
    predicted_df = pd.DataFrame(results)
    
    # save prediction results
    predicted_edges_df = pd.DataFrame([
        {
            'subject_id': edge[0]['id'],
            'subject_label': edge[0]['label'],
            'relation': edge[1]['label'],
            'object_id': edge[2]['id'],
            'object_label': edge[2]['label'],
            'notes': edge[3]['notes']
        }
        for edge in predicted_edges
    ])
    
    predicted_nodes_df = pd.DataFrame(predicted_nodes)
    
    # save to files
    predicted_edges_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_prediction_edges.csv'), index=False)
    predicted_nodes_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_prediction_nodes.csv'), index=False)
    predicted_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_pyg_predicted_df.csv'), index=False)
    
    predict_end_time = time.time()
    predict_duration = predict_end_time - predict_start_time
    formatted_predict_duration = format_duration(predict_duration)
    
    with open(pred_report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Prediction phase completed in {formatted_predict_duration} at {current_time}\n\n")
    
    # overall timing
    end_time = time.time()
    duration = end_time - start_time
    formatted_duration = format_duration(duration)
    print(f"'networkmodel_pyg.py' run finished in {formatted_duration} –where the training step took {formatted_train_duration} and the prediction {formatted_predict_duration} ({formatted_filtering_duration} of which for drug filtering and ranking).")
    
    return predicted_edges, predicted_nodes, ranked_drug_list, predicted_df