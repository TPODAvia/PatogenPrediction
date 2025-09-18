# %% [markdown]
# embed_multi_h5_3d.py — visualize multiple ProteinEmbeddings.h5
# Modes:
#  • proteins: 3D map of ALL proteins from all files (color by file; subsample per file)
#  • files:    3D map of FILE CENTROIDS (one point per file) + optional edges by cosine sim
# Also saves:  file_centroids_cosine.csv (pairwise cosine similarity between files)
#
# Dependencies: numpy, pandas, h5py, plotly (+ umap-learn optional)
# python3 embed_multi_h5_3d.py \
#   --h5_list /home/vboxuser/PatogenPrediction/h5_list.txt \
#   --mode both \
#   --method auto \
#   --sample_per_file 8000 \
#   --marker_size 1.2 \
#   --file_marker_size 8 \
#   --edges_thresh 0.92 \
#   --out_dir /home/vboxuser/PatogenPrediction/dataset/result_test/multi_h5_viz
# %%
# embed_multi_h5_3d.py — robust for few files; safe UMAP → PCA fallback

import os, argparse, numpy as np, pandas as pd, h5py
try:
    import umap
    HAS_UMAP = True
except Exception:
    HAS_UMAP = False

def ensure_3d(Y: np.ndarray) -> np.ndarray:
    """Always return (n, 3); pad with zeros if Y has <3 columns."""
    Y = np.asarray(Y)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    n, d = Y.shape
    if d == 3:
        return Y
    if d > 3:
        return Y[:, :3]
    # pad to 3
    pad = np.zeros((n, 3 - d), dtype=Y.dtype)
    return np.hstack([Y, pad])

def choose_embeddings_dataset(h5):
    cands=[]
    for k in h5.keys():
        ds=h5[k]
        if isinstance(ds,h5py.Dataset) and ds.ndim==2 and np.issubdtype(ds.dtype, np.floating):
            name=k.lower(); prefer=("emb" in name) or ("embed" in name)
            cands.append((k,prefer,ds.shape[0]*ds.shape[1]))
    if not cands: raise ValueError("No 2D float datasets in H5.")
    cands.sort(key=lambda x:(not x[1], -x[2]))
    return cands[0][0]

def choose_id_vector(h5):
    for name in ("protein_ids","protein_id","ids","names","proteins","headers"):
        if name in h5 and isinstance(h5[name],h5py.Dataset) and h5[name].ndim==1:
            arr=h5[name][:]
            try: return np.array([x.decode() if isinstance(x,(bytes,np.bytes_)) else str(x) for x in arr])
            except Exception: return np.array([str(x) for x in arr])
    return None

def sample_name_from_path(h5_path):
    parent=os.path.basename(os.path.dirname(h5_path))
    return parent.replace("preprocess_","",1) if parent.startswith("preprocess_") else os.path.splitext(os.path.basename(h5_path))[0]

def load_h5(h5_path):
    with h5py.File(h5_path,"r") as h5:
        emb_key=choose_embeddings_dataset(h5)
        X=h5[emb_key][:]; ids=choose_id_vector(h5)
    if ids is None or len(ids)!=X.shape[0]:
        ids=np.array([f"{sample_name_from_path(h5_path)}|prot_{i}" for i in range(X.shape[0])])
    return X, ids, sample_name_from_path(h5_path)

def pca_numpy(X, n_components=3, zscore=True, eps=1e-12):
    X=np.asarray(X, dtype=np.float64)
    Xc=X - X.mean(axis=0, keepdims=True)
    if zscore:
        std=Xc.std(axis=0, ddof=1, keepdims=True); std=np.where(std<eps,1.0,std); Xc=Xc/std
    U,S,Vt=np.linalg.svd(Xc, full_matrices=False)
    return (Xc @ Vt.T)[:, :n_components]

def reduce_to_3d(X, method="auto", random_state=42, n_neighbors=None):
    X=np.asarray(X, dtype=np.float64)
    n=X.shape[0]
    # Degenerate: 1–2 points → just PCA (UMAP fails)
    if n <= 2: return pca_numpy(X, 3, zscore=True)
    # UMAP path
    if method=="umap" or (method=="auto" and HAS_UMAP):
        nn = n_neighbors if n_neighbors is not None else min(15, max(2, n-1))
        if nn < 2: return pca_numpy(X, 3, zscore=True)
        try:
            reducer = umap.UMAP(n_components=3, random_state=random_state, metric="cosine", n_neighbors=nn)
            Y = reducer.fit_transform(X)
            if Y.shape[0]==0: raise ValueError("empty UMAP embedding")
            return Y
        except Exception as e:
            print(f"[warn] UMAP failed ({e}); using PCA instead.")
            return pca_numpy(X, 3, zscore=True)
    # PCA path
    return pca_numpy(X, 3, zscore=True)

def cosine_sim(A,B):
    A=np.asarray(A, dtype=np.float64); B=np.asarray(B, dtype=np.float64)
    A=A/np.maximum(np.linalg.norm(A,axis=1,keepdims=True),1e-12)
    B=B/np.maximum(np.linalg.norm(B,axis=1,keepdims=True),1e-12)
    return A @ B.T

def plot_proteins(df, out_html, marker_size=1.5):
    import plotly.express as px
    fig=px.scatter_3d(df, x="x",y="y",z="z", color="file",
                      hover_data={"protein_id":True,"file":True,"x":False,"y":False,"z":False},
                      title=f"All proteins from {df['file'].nunique()} files")
    fig.update_traces(marker=dict(size=marker_size, opacity=0.9))
    fig.update_layout(height=800, legend_title_text="File")
    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"[OK] Wrote proteins 3D HTML → {out_html}")

def plot_files(files_df, out_html, edges=None, marker_size=6.0):
    import plotly.express as px, plotly.graph_objects as go
    fig=px.scatter_3d(files_df, x="x",y="y",z="z", text="file", color="file",
                      hover_data={"file":True,"x":False,"y":False,"z":False},
                      title="File centroids (one point per H5)")
    fig.update_traces(marker=dict(size=marker_size, opacity=0.95), textposition="top center")
    if edges is not None and not edges.empty:
        xs,ys,zs=[],[],[]
        for _,r in edges.iterrows():
            a=files_df.loc[files_df["file"]==r["file_a"], ["x","y","z"]].values[0]
            b=files_df.loc[files_df["file"]==r["file_b"], ["x","y","z"]].values[0]
            xs += [a[0],b[0],None]; ys += [a[1],b[1],None]; zs += [a[2],b[2],None]
        fig.add_trace(go.Scatter3d(x=xs,y=ys,z=zs, mode="lines", line=dict(width=2),
                                   name=f"Edges (cos≥{edges.attrs.get('thresh','…')})",
                                   hoverinfo="none", showlegend=True))
    fig.update_layout(height=800, legend_title_text="File")
    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"[OK] Wrote file-centroid 3D HTML → {out_html}")

def main():
    ap=argparse.ArgumentParser(description="3D correlation view across multiple ProteinEmbeddings.h5")
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--h5", nargs="+")
    g.add_argument("--h5_list")
    ap.add_argument("--mode", choices=["proteins","files","both"], default="both")
    ap.add_argument("--method", choices=["auto","umap","pca"], default="auto")
    ap.add_argument("--sample_per_file", type=int, default=10000)
    ap.add_argument("--marker_size", type=float, default=1.2)
    ap.add_argument("--file_marker_size", type=float, default=8.0)
    ap.add_argument("--edges_thresh", type=float, default=0.9)
    ap.add_argument("--out_dir", required=True)
    args=ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    paths=[]
    if args.h5_list:
        with open(args.h5_list) as f:
            paths=[ln.strip() for ln in f if ln.strip()]
    else:
        paths=args.h5
    paths=[p for p in paths if os.path.isfile(p)]
    if not paths: raise SystemExit("No valid H5 files provided.")

    all_points=[]; centroids=[]; file_names=[]
    for p in paths:
        X, ids, label = load_h5(p)
        file_names.append(label)
        centroids.append(np.mean(X, axis=0, dtype=np.float64))
        if args.mode in ("proteins","both"):
            if args.sample_per_file and args.sample_per_file>0 and X.shape[0]>args.sample_per_file:
                rng=np.random.default_rng(42)
                idx=np.sort(rng.choice(X.shape[0], size=args.sample_per_file, replace=False))
                Xp=X[idx]; idsp=ids[idx]
            else:
                Xp=X; idsp=ids
            all_points.append((Xp, idsp, np.full(Xp.shape[0], label, dtype=object)))

    if args.mode in ("proteins","both"):
        Xcat=np.vstack([t[0] for t in all_points])
        labels=np.concatenate([t[2] for t in all_points])
        prot_ids=np.concatenate([t[1] for t in all_points])
        Y3=reduce_to_3d(Xcat, method=args.method) 
        Y3 = ensure_3d(Y3)
        dfp=pd.DataFrame({"x":Y3[:,0],"y":Y3[:,1],"z":Y3[:,2],"file":labels,"protein_id":prot_ids})
        plot_proteins(dfp, os.path.join(args.out_dir,"proteins_3D.html"), marker_size=args.marker_size)

    if args.mode in ("files","both"):
        C=np.vstack(centroids)
        S=cosine_sim(C,C)
        pd.DataFrame(S, index=file_names, columns=file_names).to_csv(os.path.join(args.out_dir,"file_centroids_cosine.csv"))
        n_files=C.shape[0]
        # For few files, force PCA; else safe-UMAP with bounded n_neighbors
        method_files = "pca" if (args.method=="pca" or n_files<5 or not HAS_UMAP and args.method!="umap") else args.method
        nn = None
        if method_files in ("umap","auto") and HAS_UMAP:
            method_files="umap"
            nn=min(15, max(2, n_files-1))
        Yf=reduce_to_3d(C, method=method_files, n_neighbors=nn)
        Yf = ensure_3d(Yf)
        dff=pd.DataFrame({"x":Yf[:,0],"y":Yf[:,1],"z":Yf[:,2],"file":file_names})
        edges=[]
        if args.edges_thresh and args.edges_thresh>0:
            for i in range(n_files):
                for j in range(i+1, n_files):
                    if S[i,j] >= args.edges_thresh:
                        edges.append({"file_a":file_names[i], "file_b":file_names[j], "sim":S[i,j]})
        edges_df=pd.DataFrame(edges)
        if not edges_df.empty: edges_df.attrs["thresh"]=args.edges_thresh
        plot_files(dff, os.path.join(args.out_dir,"files_3D.html"),
                   edges=edges_df if not edges_df.empty else None,
                   marker_size=args.file_marker_size)

if __name__=="__main__":
    main()