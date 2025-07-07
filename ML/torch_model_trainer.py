import torch
import torch.nn.functional as F
from torch_geometric.nn import GINConv, GINEConv, GATConv, global_add_pool
from ML.torch_data_processing import read_targets, load_data_from_sdf, create_dataloader
import os

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
        self.node_encoder = torch.nn.Linear(input_dim, hidden_dim)

        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            mlp = torch.nn.Sequential(
                torch.nn.Linear(hidden_dim, hidden_dim),
                torch.nn.ReLU(),
                torch.nn.Linear(hidden_dim, hidden_dim)
            )
            conv = GINEConv(mlp, edge_dim=edge_dim)
            self.convs.append(conv)

        self.readout = global_add_pool
        self.output = torch.nn.Linear(hidden_dim, 1)

    def forward(self, x, edge_index, edge_attr, batch):
        x = self.node_encoder(x)
        for conv in self.convs:
            x = conv(x, edge_index, edge_attr)
            x = F.relu(x)
        x = self.readout(x, batch)
        return self.output(x).squeeze()


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


def train_and_save_model(
    sdf_dir,
    target_file,
    modelo_nombre,
    epochs,
    save_path,
    batch_size=32,
    lr=0.001
):
    target_dict = read_targets(target_file)
    targetname = os.path.splitext(os.path.basename(target_file))[0]
    data_list = load_data_from_sdf(sdf_dir, target_dict)
    dataloader = create_dataloader(data_list, batch_size=batch_size)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    input_dim = data_list[0].x.shape[1]
    edge_dim = data_list[0].edge_attr.shape[1]
    
    model = create_model(modelo_nombre, input_dim, edge_dim)

    train(model, dataloader, device, epochs=epochs, lr=lr)

    checkpoint = {
        'model_state_dict': model.state_dict(),
        'model_type': modelo_nombre,
        'input_dim': input_dim,
        'edge_dim': edge_dim,
        'epochs_trained': epochs,
        'target_name': targetname,
    }

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    torch.save(checkpoint, save_path)

    return save_path  # opcional, para confirmar ruta

