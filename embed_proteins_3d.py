# %% [markdown]
# embed_proteins_3d.py — small markers + legend + gene location hover
# - Pulls contig/start/end/strand/aa_len from PredictedProteins.faa (same folder as H5)
# - Color by contig (or strand/length), symbol by strand
# - No SciPy/sklearn required (UMAP optional; fallback to NumPy PCA)

'''
python3 embed_multi_h5_3d.py   --h5_list /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/h5_list.txt   --mode both   --method auto   --sample_per_file 8000   --marker_size 1.2   --file_marker_size 8   --edges_thresh 0.92   --out_dir /home/rover2/HW1_popgen/kursov/PatogenPrediction/multi_h5_viz
'''

# %%
import os
import argparse
import numpy as np
import pandas as pd
import h5py

# Try UMAP (optional)
try:
    import umap
    HAS_UMAP = True
except Exception:
    HAS_UMAP = False

# ---------- helpers ----------
def choose_embeddings_dataset(h5: h5py.File) -> str:
    cands = []
    for k in h5.keys():
        ds = h5[k]
        if isinstance(ds, h5py.Dataset) and ds.ndim == 2 and np.issubdtype(ds.dtype, np.floating):
            name = k.lower()
            prefer = ("emb" in name) or ("embed" in name)
            cands.append((k, prefer, ds.shape[0]*ds.shape[1]))
    if not cands:
        raise ValueError("No 2D float datasets found in H5.")
    cands.sort(key=lambda x: (not x[1], -x[2]))
    return cands[0][0]

def choose_id_vector(h5: h5py.File):
    for name in ("protein_ids","protein_id","ids","names","proteins","headers"):
        if name in h5 and isinstance(h5[name], h5py.Dataset) and h5[name].ndim == 1:
            arr = h5[name][:]
            try:
                return np.array([x.decode() if isinstance(x,(bytes,np.bytes_)) else str(x) for x in arr])
            except Exception:
                return np.array([str(x) for x in arr])
    return None

def parse_prodigal_faa(faa_path: str):
    """
    Parse PredictedProteins.faa from Prodigal.
    Header pattern: '>contig # start # end # strand # ID=...;...' (with spaces around '#')
    Returns dict of arrays: contig, start, end, strand, aa_len, header
    """
    contigs, starts, ends, strands, aalen, headers = [], [], [], [], [], []
    if not os.path.isfile(faa_path):
        return None
    cur_len = 0
    cur_header = None
    def flush():
        if cur_header is not None:
            headers.append(cur_header)
            aalen.append(cur_len)

    with open(faa_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                cur_len = 0
                h = line[1:].strip()
                headers.append(h)  # temp; will overwrite in flush
                headers.pop()      # keep flow simple
                # parse contig/start/end/strand
                parts = h.split(" # ")
                contig = parts[0].split()[0] if parts else None
                try:
                    start = int(parts[1].split()[0]) if len(parts) > 1 else None
                    end   = int(parts[2].split()[0]) if len(parts) > 2 else None
                    st    = parts[3].split()[0] if len(parts) > 3 else None
                except Exception:
                    start = end = st = None
                strand = { "1":"+","-1":"-","+":"+", "-":"-" }.get(str(st), None)
                contigs.append(contig); starts.append(start); ends.append(end); strands.append(strand)
                cur_header = h
            else:
                cur_len += len(line.strip())
        flush()

    return {
        "contig": np.array(contigs, dtype=object),
        "start":  np.array(starts, dtype=object),
        "end":    np.array(ends, dtype=object),
        "strand": np.array(strands, dtype=object),
        "aa_len": np.array(aalen, dtype=object),
        "header": np.array(headers, dtype=object),
    }

def pca_numpy(X: np.ndarray, n_components: int = 3, zscore: bool = True, eps: float = 1e-12) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    Xc = X - X.mean(axis=0, keepdims=True)
    if zscore:
        std = Xc.std(axis=0, ddof=1, keepdims=True)
        std = np.where(std < eps, 1.0, std)
        Xc = Xc / std
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    scores = Xc @ Vt.T
    return scores[:, :n_components]

def reduce_to_3d(X: np.ndarray, method: str = "auto", random_state: int = 42) -> np.ndarray:
    if method == "umap" or (method == "auto" and HAS_UMAP):
        reducer = umap.UMAP(n_components=3, random_state=random_state, metric="cosine")
        return reducer.fit_transform(X)
    return pca_numpy(X, n_components=3, zscore=True)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="3D interactive embedding of PF2 proteins with gene-location hover")
    ap.add_argument("--h5", required=True, help="Path to ProteinEmbeddings.h5")
    ap.add_argument("--out_html", default=None, help="Output HTML (default: alongside H5)")
    ap.add_argument("--method", choices=["auto","umap","pca"], default="auto")
    ap.add_argument("--sample", type=int, default=20000, help="Subsample N proteins for speed (0 = all)")
    ap.add_argument("--marker_size", type=float, default=2.0, help="Marker size (default: 2.0, smaller = finer)")
    ap.add_argument("--color_by", choices=["contig","strand","length","none"], default="contig",
                    help="Legend/coloring field")
    args = ap.parse_args()

    h5_path = args.h5
    if args.out_html is None:
        stem = os.path.splitext(os.path.basename(h5_path))[0]
        args.out_html = os.path.join(os.path.dirname(h5_path), f"{stem}_3D.html")

    # Load embeddings + ids
    with h5py.File(h5_path, "r") as h5:
        emb_key = choose_embeddings_dataset(h5)
        X = h5[emb_key][:]
        ids = choose_id_vector(h5)
        if ids is None or len(ids) != X.shape[0]:
            ids = np.array([f"prot_{i}" for i in range(X.shape[0])])

    # Pull gene metadata from sibling PredictedProteins.faa
    preprocess_dir = os.path.dirname(h5_path)
    faa_path = os.path.join(preprocess_dir, "PredictedProteins.faa")
    meta = parse_prodigal_faa(faa_path)
    if (meta is None) or (len(meta["contig"]) != X.shape[0]):
        print("[warn] Could not align FAA metadata to embeddings — hover will only show protein_id.")
        meta = {
            "contig": np.array([None]*X.shape[0], dtype=object),
            "start":  np.array([None]*X.shape[0], dtype=object),
            "end":    np.array([None]*X.shape[0], dtype=object),
            "strand": np.array([None]*X.shape[0], dtype=object),
            "aa_len": np.array([None]*X.shape[0], dtype=object),
            "header": np.array([None]*X.shape[0], dtype=object),
        }

    # Subsample (keep metadata in sync)
    N = X.shape[0]
    if args.sample and args.sample > 0 and N > args.sample:
        rng = np.random.default_rng(42)
        idx = np.sort(rng.choice(N, size=args.sample, replace=False))
        X = X[idx]; ids = ids[idx]
        for k in meta:
            meta[k] = meta[k][idx]
        print(f"[info] Subsampled {args.sample} / {N} proteins.")
    else:
        idx = np.arange(N)

    # Compute 3D coords
    method = "pca" if (args.method == "auto" and not HAS_UMAP) else args.method
    Y3 = reduce_to_3d(X, method=method)

    # Prepare plot dataframe
    df = pd.DataFrame({
        "x": Y3[:,0], "y": Y3[:,1], "z": Y3[:,2],
        "protein_id": ids,
        "contig": meta["contig"],
        "start": meta["start"],
        "end": meta["end"],
        "strand": meta["strand"],
        "aa_len": meta["aa_len"],
    })
    # optional length bins for nicer legend
    df["length_bin"] = pd.cut(pd.to_numeric(df["aa_len"], errors="coerce"),
                              bins=[0,100,200,400,800, 20000],
                              labels=["≤100","100–200","200–400","400–800",">800"])

    # Choose color mapping
    color_col = None
    if args.color_by == "contig":
        color_col = "contig"
    elif args.color_by == "strand":
        color_col = "strand"
    elif args.color_by == "length":
        color_col = "length_bin"

    # Build plot
    import plotly.express as px
    fig = px.scatter_3d(
        df, x="x", y="y", z="z",
        color=color_col if color_col else None,
        symbol="strand" if "strand" in df.columns else None,
        hover_data={
            "protein_id": True,
            "contig": True,
            "start": True,
            "end": True,
            "strand": True,
            "aa_len": True,
            "x": False, "y": False, "z": False
        },
        title=f"3D protein embedding ({method.upper()}): {os.path.basename(h5_path)}"
    )
    # smaller dots + legend shown
    fig.update_traces(marker=dict(size=args.marker_size, opacity=0.9))
    fig.update_layout(
        height=800,
        legend_title_text=f"Color: {args.color_by}",
        legend=dict(itemsizing="constant")
    )

    fig.write_html(args.out_html, include_plotlyjs="cdn")
    print(f"[OK] Wrote interactive 3D HTML → {args.out_html}")

if __name__ == "__main__":
    main()