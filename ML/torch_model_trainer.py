import torch
import torch.nn.functional as F
from torch_geometric.nn import GINConv, GINEConv, GATConv, global_add_pool
from torch_geometric.data import DataLoader

from torch import nn
from torch_geometric.nn import Sequential

# ----------------------
# Modelos GNN
# ----------------------

class GINNet(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim=64, num_layers=3):
        super().__init__()
        # Construimos capas GIN con MLP interno
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            mlp = torch.nn.Sequential(
                torch.nn.Linear(input_dim if i == 0 else hidden_dim, hidden_dim),
                torch.nn.ReLU(),
                torch.nn.Linear(hidden_dim, hidden_dim)
            )
            conv = GINConv(mlp)
            self.convs.append(conv)
        self.lin = torch.nn.Linear(hidden_dim, 1)  # Regresión

    def forward(self, x, edge_index, edge_attr, batch):
        for conv in self.convs:
            x = conv(x, edge_index)
            x = F.relu(x)
        x = global_add_pool(x, batch)  # Pooling global por grafo
        out = self.lin(x)
        return out.squeeze()


class GINENet(torch.nn.Module):
    def __init__(self, input_dim, edge_dim=1, hidden_dim=64, num_layers=3):
        super().__init__()
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            nn_edge = torch.nn.Sequential(
                torch.nn.Linear(edge_dim, hidden_dim),
                torch.nn.ReLU(),
                torch.nn.Linear(hidden_dim, hidden_dim)
            )
            conv = GINEConv(nn_edge)
            self.convs.append(conv)
        self.lin = torch.nn.Linear(hidden_dim, 1)

    def forward(self, x, edge_index, edge_attr, batch):
        for conv in self.convs:
            x = conv(x, edge_index, edge_attr)
            x = F.relu(x)
        x = global_add_pool(x, batch)
        out = self.lin(x)
        return out.squeeze()


class GATNet(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim=64, num_layers=3, heads=4):
        super().__init__()
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            in_channels = input_dim if i == 0 else hidden_dim * heads
            conv = GATConv(in_channels, hidden_dim, heads=heads, concat=True)
            self.convs.append(conv)
        self.lin = torch.nn.Linear(hidden_dim * heads, 1)

    def forward(self, x, edge_index, edge_attr, batch):
        for conv in self.convs:
            x = conv(x, edge_index)
            x = F.elu(x)
        x = global_add_pool(x, batch)
        out = self.lin(x)
        return out.squeeze()

# ----------------------
# Función para crear modelo según elección
# ----------------------

def create_model(model_name, input_dim, edge_dim=1):
    if model_name == "GIN":
        return GINNet(input_dim=input_dim)
    elif model_name == "GINE":
        return GINENet(input_dim=input_dim, edge_dim=edge_dim)
    elif model_name == "GAT":
        return GATNet(input_dim=input_dim)
    else:
        raise ValueError(f"Modelo desconocido: {model_name}")

# ----------------------
# Función para entrenar modelo
# ----------------------

def train(model, dataloader, device, epochs=20, lr=0.001):
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = torch.nn.MSELoss()  # Regresión

    model.train()
    for epoch in range(1, epochs + 1):
        total_loss = 0
        for batch in dataloader:
            batch = batch.to(device)
            optimizer.zero_grad()
            # Los modelos reciben (x, edge_index, edge_attr, batch)
            out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
            loss = criterion(out, batch.y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * batch.num_graphs
        print(f"Epoch {epoch:03d} | Loss: {total_loss / len(dataloader.dataset):.6f}")

