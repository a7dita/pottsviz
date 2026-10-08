#!/usr/bin/env python3
"""Cache the notebook's display-only spectral permutations (requires NumPy).

Matches simcode_pilot.analysis.topology_analysis.plot_adjacency_atlas.
The CSVs, original memory IDs and simulation matrices are never rewritten.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[1]
directory = root / "src/lib/topologies"
manifest = json.loads((directory / "manifest.json").read_text())
orders = {}
for topology in manifest["topologies"]:
    content = (directory / topology["file"]).read_bytes()
    matrix = np.loadtxt(directory / topology["file"], delimiter=",")
    frontal = posterior = np.arange(len(matrix))
    if topology["degree"] == 7:
        u, _, vt = np.linalg.svd(matrix.astype(float), full_matrices=False)
        frontal = np.argsort(u[:, 1])
        posterior = np.argsort(vt[1])
    orders[topology["id"]] = {
        "source_blob": hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest(),
        "frontal": frontal.tolist(),
        "posterior": posterior.tolist(),
    }
(root / "src/lib/topology-display-order.json").write_text(json.dumps(orders, indent=2) + "\n")
print(f"Saved notebook display orders for {len(orders)} topologies.")
