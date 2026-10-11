#!/usr/bin/env python3
"""Generate 17 ordered structural association conditions and fresh metrics.

All CSVs are binary frontal-row × posterior-column matrices. Generator labels
are retained: no row/column relabeling, shuffling, or spectral display ordering.
Groups are structural, with no assumed semantic or spatial interpretation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.csgraph import connected_components

M = 49
K = 7
DEFAULT_SEED = 20261009
METRIC_KEYS = ("redundancy_log1p", "q_star", "mixing_gap")


def component_labels(graph: np.ndarray) -> np.ndarray:
    zero = csr_matrix((M, M))
    adjacency = bmat([[zero, csr_matrix(graph)], [csr_matrix(graph.T), zero]], format="csr")
    return connected_components(adjacency, directed=False, return_labels=True)[1]


def is_connected(graph: np.ndarray) -> bool:
    return len(np.unique(component_labels(graph))) == 1


def validate_graph(graph: np.ndarray, require_connected: bool = True, degree: int = K) -> None:
    if graph.shape != (M, M) or not np.isin(graph, [0, 1]).all():
        raise ValueError("expected a binary 49 × 49 matrix")
    if not np.all(graph.sum(axis=1) == degree) or not np.all(graph.sum(axis=0) == degree):
        raise ValueError(f"every memory must have exactly {degree} partners")
    if require_connected and not is_connected(graph):
        raise ValueError("topology is disconnected")


def block_graph() -> np.ndarray:
    graph = np.zeros((M, M), dtype=np.uint8)
    for group in range(K):
        block = slice(group * K, (group + 1) * K)
        graph[block, block] = 1
    return graph


def modular_graph(m: int, seed: int = DEFAULT_SEED) -> np.ndarray:
    """Exactly m partners in the corresponding group, on both sides."""
    if not 0 <= m <= K:
        raise ValueError("m must be between zero and seven")
    rng = np.random.default_rng(seed)
    group = np.arange(M) // K
    same = group[:, None] == group[None, :]
    graph = np.zeros((M, M), dtype=np.uint8)
    for within in [True] * m + [False] * (K - m):
        allowed = (same if within else ~same) & (graph == 0)
        cost = rng.random((M, M))
        cost[~allowed] = 1e6
        rows, cols = linear_sum_assignment(cost)
        if not allowed[rows, cols].all():
            raise RuntimeError("unable to construct a disjoint allowed perfect matching")
        graph[rows, cols] = 1
    validate_graph(graph, require_connected=m < K)
    return graph


def ring_neighborhood() -> np.ndarray:
    """Seven local targets: cyclic offsets 0, ±1, ±2, ±3."""
    indices = np.arange(M)
    offsets = (indices[None, :] - indices[:, None]) % M
    return np.isin(offsets, [0, 1, 2, 3, M - 3, M - 2, M - 1])


def ring_graph(r: int, seed: int = DEFAULT_SEED) -> np.ndarray:
    """Exactly r local and 7-r outside partners per memory, on both sides.

    r counts local connections, not jump distance or pairwise shared targets.
    Matching choices are random within the allowed regions, as for modular_graph.
    """
    if not 1 <= r <= K:
        raise ValueError("r must be between one and seven")
    rng = np.random.default_rng(seed)
    local = ring_neighborhood()
    graph = np.zeros((M, M), dtype=np.uint8)
    for within in [True] * r + [False] * (K - r):
        allowed = (local if within else ~local) & (graph == 0)
        cost = rng.random((M, M))
        cost[~allowed] = 1e6
        rows, cols = linear_sum_assignment(cost)
        if not allowed[rows, cols].all():
            raise RuntimeError("unable to construct a disjoint allowed perfect matching")
        graph[rows, cols] = 1
    validate_graph(graph)
    if not all(np.all((graph * local).sum(axis=axis) == r) for axis in (0, 1)):
        raise RuntimeError("every memory must have exactly r local partners")
    return graph


def random_regular(seed: int) -> np.ndarray:
    """Degree-preserving randomization of edges, without relabeling memories."""
    rng = np.random.default_rng(seed)
    graph = np.zeros((M, M), dtype=np.uint8)
    for row in range(M):
        graph[row, (row + np.arange(K)) % M] = 1
    accepted = 0
    while accepted < 40000:
        f1, f2 = rng.choice(M, size=2, replace=False)
        p1 = int(rng.choice(np.flatnonzero(graph[f1])))
        p2 = int(rng.choice(np.flatnonzero(graph[f2])))
        if p1 == p2 or graph[f1, p2] or graph[f2, p1]:
            continue
        graph[f1, p1] = graph[f2, p2] = 0
        graph[f1, p2] = graph[f2, p1] = 1
        accepted += 1
    validate_graph(graph)
    return graph


def redundancy_count(graph: np.ndarray) -> int:
    shared = graph.astype(np.int64) @ graph.astype(np.int64).T
    upper = shared[np.triu_indices(M, k=1)]
    return int(np.sum(upper * (upper - 1) // 2))


def mixing_gap(graph: np.ndarray) -> float:
    singular_values = np.linalg.svd(graph.astype(float) / graph.sum(axis=1)[0], compute_uv=False)
    return float(np.clip(1.0 - singular_values[1] ** 2, 0.0, 1.0))


def barber_modularity(graph: np.ndarray, labels_f: np.ndarray, labels_p: np.ndarray) -> float:
    degree = int(graph.sum(axis=1)[0])
    expected = degree / M
    same = labels_f[:, None] == labels_p[None, :]
    return float(np.sum((graph - expected) * same) / (M * degree))


def _coordinate_ascent_modularity(
    graph: np.ndarray, labels_f: np.ndarray, labels_p: np.ndarray, rng: np.random.Generator
) -> tuple[float, np.ndarray, np.ndarray]:
    labels_f = labels_f.copy()
    labels_p = labels_p.copy()
    degree = int(graph.sum(axis=1)[0])
    edge_total = float(M * degree)
    for _ in range(40):
        moved = False
        for side in rng.permutation(2):
            order = rng.permutation(M)
            for node in order:
                if side == 0:
                    current = int(labels_f[node])
                    candidates = np.unique(labels_p)
                    d_other = np.bincount(labels_p, minlength=2 * M) * degree
                    edge_by_label = np.bincount(labels_p, weights=graph[node], minlength=2 * M)
                else:
                    current = int(labels_p[node])
                    candidates = np.unique(labels_f)
                    d_other = np.bincount(labels_f, minlength=2 * M) * degree
                    edge_by_label = np.bincount(labels_f, weights=graph[:, node], minlength=2 * M)
                base = edge_by_label[current] - degree * d_other[current] / edge_total
                gains = edge_by_label[candidates] - degree * d_other[candidates] / edge_total - base
                best_pos = int(np.argmax(gains))
                if gains[best_pos] > 1e-12:
                    target = int(candidates[best_pos])
                    if side == 0:
                        labels_f[node] = target
                    else:
                        labels_p[node] = target
                    moved = True
        if not moved:
            break
    return barber_modularity(graph, labels_f, labels_p), labels_f, labels_p


def q_star(graph: np.ndarray, seed: int = 0) -> float:
    """Fixed-restart greedy estimate of maximum Barber bipartite modularity."""
    rng = np.random.default_rng(seed)
    # Component partitions are exact for unions of complete bipartite blocks.
    labels = component_labels(graph)
    best = barber_modularity(graph, labels[:M], labels[M:])
    if all(np.all(graph[np.ix_(np.flatnonzero(labels[:M] == c),
                                 np.flatnonzero(labels[M:] == c))])
           for c in np.unique(labels)):
        return best
    # Include the structural partition as a reproducible feasible starting point.
    planted = np.arange(M) // K
    value, _, _ = _coordinate_ascent_modularity(graph, planted, planted, rng)
    best = max(best, value)
    for communities in (2, 3, 4, 5, 7, 10, 14):
        for _ in range(3):
            labels_f = rng.integers(0, communities, size=M)
            labels_p = rng.integers(0, communities, size=M)
            value, _, _ = _coordinate_ascent_modularity(graph, labels_f, labels_p, rng)
            best = max(best, value)
    return float(best)


def graph_metrics(graph: np.ndarray) -> dict[str, float | int]:
    n4 = redundancy_count(graph)
    return {
        "n4": n4,
        "redundancy_log1p": float(math.log1p(n4)),
        "q_star": q_star(graph, seed=1729),
        "mixing_gap": mixing_gap(graph),
    }


def save_matrix(path: Path, graph: np.ndarray) -> None:
    path.write_text("\n".join(",".join(str(int(v)) for v in row) for row in graph) + "\n", encoding="utf-8")


def design_graphs(seed: int = DEFAULT_SEED) -> list[tuple[str, str, int | None, np.ndarray]]:
    return [("original", "original", None, np.eye(M, dtype=np.uint8)),
            ("random", "random", None, random_regular(seed))] + [
        (f"m{m}", "modular", m, modular_graph(m, seed + 100 + m)) for m in range(K)
    ] + [("m7", "modular", K, block_graph())] + [
        (f"r{r}", "ring-local", r, ring_graph(r, seed + 200 + r)) for r in range(1, K + 1)
    ]


def build_design(output: Path, seed: int = DEFAULT_SEED) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    previous_manifest = output / "manifest.json"
    previous_files = set()
    if previous_manifest.exists():
        previous_files = {item["file"] for item in json.loads(previous_manifest.read_text())["topologies"]}
    entries = []
    for name, family, parameter, graph in design_graphs(seed):
        degree = 1 if name == "original" else K
        validate_graph(graph, require_connected=name not in ("original", "m7"), degree=degree)
        filename = f"{name}.csv"
        save_matrix(output / filename, graph)
        metric = graph_metrics(np.loadtxt(output / filename, delimiter=",", dtype=np.uint8))
        entries.append({
            "id": name, "name": name, "family": family, "parameter": parameter,
            "role": "modular_endpoint" if name == "m7" else {
                "original": "one_to_one_reference", "random": "random_regular_reference",
            }.get(family, family),
            "file": filename, "memories_per_half": M, "degree": degree,
            "edge_weight_after_loading": 1.0 / degree,
            "components": int(len(np.unique(component_labels(graph)))),
            "metrics": metric,
        })
    manifest = {
        "schema_version": 4, "generation_seed": seed, "label_order": "generator_order",
        "matrix_axes": {"rows": "frontal", "columns": "posterior"},
        "selection": {"description": "original, random, m0–m7, and ring-local r1–r7.",
                      "coordinates": list(METRIC_KEYS)},
        "metric_method": {"q_star": "Barber modularity: component and structural starts, plus 21 random starts; approximate except complete components",
                          "q_star_seed": 1729, "mixing_gap": "1 - sigma2(B / degree)^2"},
        "families": {"modular": {"parameter": "m", "min": 0, "max": 7, "endpoint": "m7"},
                     "ring-local": {"parameter": "r", "min": 1, "max": 7, "endpoint": "r7",
                                    "local_offsets": [0, -1, 1, -2, 2, -3, 3],
                                    "seed_offset": 200}},
        "topologies": entries,
    }
    # Replace the previous design while retaining unrelated files and documentation.
    retained = {item["file"] for item in entries}
    for filename in previous_files - retained:
        old_csv = output / filename
        if old_csv.parent == output and old_csv.suffix == ".csv" and old_csv.exists():
            old_csv.unlink()
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    with (output / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=["topology", "family", "parameter", "degree", "components", "n4", *METRIC_KEYS])
        writer.writeheader()
        for item in entries:
            writer.writerow({"topology": item["id"], "family": item["family"], "parameter": item["parameter"],
                             "degree": item["degree"], "components": item["components"], **item["metrics"]})
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    repository = Path(__file__).resolve().parents[1]
    default_output = repository / "src/lib/topologies" if (repository / "src/lib").is_dir() else repository / "topologies"
    parser.add_argument("--output", type=Path, default=default_output)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    manifest = build_design(args.output, args.seed)
    for item in manifest["topologies"]:
        metric = item["metrics"]
        print(f"{item['id']:8} k={item['degree']} N4={metric['n4']:4} "
              f"R4={metric['redundancy_log1p']:.4f} Q*={metric['q_star']:.4f} gamma={metric['mixing_gap']:.4f}")


if __name__ == "__main__":
    main()
