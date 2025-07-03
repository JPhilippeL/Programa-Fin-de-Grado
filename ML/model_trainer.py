# model_trainer.py

from tensorflow import keras
import tensorflow as tf
from molgraph.chemistry import MolecularGraphEncoder
from molgraph import layers
import numpy as np

def build_model(model_type: str, input_spec, hidden_units=32):
    """
    Construye un modelo GNN basado en molgraph.

    Args:
      model_type: "mpnn", "gin", o "gat"
      input_spec: GraphTensorSpec del encoder
      hidden_units: unidades en capas ocultas

    Devuelve:
      modelo compilado (Keras Model)
    """
    if model_type not in {"mpnn", "gin", "gat"}:
        raise ValueError("model_type debe ser 'mpnn', 'gin' o 'gat'")

    inputs = keras.layers.Input(type_spec=input_spec)

    if model_type == "mpnn":
        x = layers.MPNN(units=hidden_units)(inputs)
    elif model_type == "gin":
        x = layers.GINConv(units=hidden_units)(inputs)
        x = layers.GINConv(units=hidden_units)(x)
    elif model_type == "gat":
        x = layers.GATv2Conv(units=hidden_units)(inputs)
        x = layers.GATv2Conv(units=hidden_units)(x)

    x = layers.Readout()(x)
    outputs = keras.layers.Dense(units=1)(x)

    model = keras.Model(inputs, outputs)
    model.compile(optimizer="adam", loss="mse")
    return model

def entrenar_modelo(sdf_dir, target_file, model_type="mpnn", epochs=20, batch_size=8):
    """
    Carga moléculas, entrena un modelo GNN y devuelve el modelo entrenado.

    Args:
      sdf_dir: carpeta con archivos .sdf
      target_file: archivo con nombres y valores objetivo
      model_type: "mpnn", "gin" o "gat"
      epochs, batch_size: hiperparámetros

    Devuelve:
      modelo entrenado (Keras Model)
    """
    # Encoder
    encoder = MolecularGraphEncoder(atom_features="basic", bond_features="basic")
    # Cargar datos
    molgraphs, names = [], []
    import os
    from rdkit import Chem

    for fname in sorted(os.listdir(sdf_dir)):
        if fname.endswith(".sdf"):
            mol = Chem.MolFromMolFile(os.path.join(sdf_dir, fname), sanitize=True)
            if mol:
                molgraphs.append(encoder(mol))
                names.append(os.path.splitext(fname)[0])

    # Cargar objetivos
    target_map = {}
    with open(target_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts)==2:
                target_map[parts[0]] = float(parts[1])
    y = np.array([target_map.get(n, 0.0) for n in names])

    # Crear dataset
    ds = tf.data.Dataset.from_tensor_slices((molgraphs, y))
    ds = ds.shuffle(len(molgraphs)).batch(batch_size)

    # Construir modelo
    model = build_model(model_type, molgraphs[0].spec)

    # Entrenar
    model.fit(ds, epochs=epochs)
    return model
