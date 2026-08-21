#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run_crossNN.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# DNA methylation-based CNS tumour classification with crossNN
# (Yuan et al., Nature Cancer 2025; https://gitlab.com/euskirchen-lab/crossNN)
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: envs/crossNN.yaml
# Usage:
"""
python3 scripts/Run_crossNN.py \
        -i {input.adata} \
        -m {params.model_dir} \
        -w {params.weights} \
        -b {params.threshold} \
        -s {params.score_cutoff} \
        -n {params.min_features} \
        -t {threads} \
        -o_pred {output.predictions} \
        -o_scores {output.scores} \
        -o_coverage {output.coverage}
"""
#
# TODO:
# 1) Validate score scaling on cases with known integrated diagnosis
#
# History:
#  18-08-2026: File creation
#  18-08-2026: Load bundled .pkl model (Capper_et_al_NN.pkl)
#  18-08-2026: Recursive harvest of nested pickle; rescale-to-full logits
#  18-08-2026: Reach probe IDs held in a DataFrame index/column
#  20-08-2026: FIX - inputs are ternary (+1/-1/0) at threshold 0.623, not beta
#              values. Earlier versions passed raw betas, which encoded every
#              unmethylated site as the training-time mask value and produced
#              wrong class rankings. Rescaling for missing probes removed:
#              masked features are 0 by design.
#  20-08-2026: FIX - logits are z-scored across classes before the softmax, and
#              the threshold is 0.6. Both verified against upstream
#              NN_model.py (NN_classifier.predict) and training.py rather than
#              inferred from the paper.
#
# Model layout (Capper_et_al_NN.pkl), confirmed by --inspect:
#   [0] dict {'layer_out.weight': tensor (91, 366263)}   bias-free linear layer
#   [1] sklearn LabelEncoder, .classes_ = 91 methylation classes
#   [2] DataFrame (366263, 1) holding the CpG feature space
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Run_crossNN.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Classify CNS tumours from methylation beta values using crossNN  ')
    parser.add_argument('-i', help='path to preprocessed methylation (.h5ad)',
                        dest='adata',
                        type=str)
    parser.add_argument('-m', help='path to cloned crossNN repository',
                        dest='model_dir',
                        type=str)
    parser.add_argument('-w', help='absolute path to the model pickle',
                        dest='weights',
                        type=str)
    parser.add_argument('-c', help='optional class label file (one per line); '
                                   'only needed if absent from the pickle',
                        dest='classes',
                        type=str, default=None)
    parser.add_argument('-f', help='optional feature (CpG) file (one per line); '
                                   'only needed if absent from the pickle',
                        dest='features',
                        type=str, default=None)
    parser.add_argument('-l', help='adata layer holding betas (default: X)',
                        dest='layer',
                        type=str, default=None)
    parser.add_argument('-b', help='beta binarization threshold (default 0.6, '
                                   'per training.py; the paper states 0.623)',
                        dest='threshold',
                        type=float, default=0.6)
    parser.add_argument('-s', help='confidence score cutoff (array: 0.4)',
                        dest='score_cutoff',
                        type=float, default=0.4)
    parser.add_argument('-n', help='minimum usable CpGs per sample',
                        dest='min_features',
                        type=int, default=5000)
    parser.add_argument('-t', help='number of threads',
                        dest='threads',
                        type=int, default=1)
    parser.add_argument('--inspect', help='dump pickle structure and exit',
                        dest='inspect',
                        action='store_true')
    parser.add_argument('-o_pred', help='path to per-sample predictions (.tsv)',
                        dest='predictions',
                        type=str)
    parser.add_argument('-o_scores', help='path to full score matrix (.tsv)',
                        dest='scores',
                        type=str)
    parser.add_argument('-o_coverage', help='path to feature coverage report (.tsv)',
                        dest='coverage',
                        type=str)
    args = parser.parse_args()
    return args
args = parse_args()
"""
args.adata = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/methylation_data.h5ad'
args.model_dir = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/crossNN/models'
args.weights = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/crossNN/models/models/Capper_et_al_NN.pkl'
args.classes = None
args.features = None
args.layer = None
args.score_cutoff = 0.4
args.min_features = 5000
args.threads = 2
args.predictions = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/crossNN/crossNN_predictions.tsv'
args.scores = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/crossNN/crossNN_scores.tsv'
args.coverage = '/trinity/home/r115502/SSLOWGRADE/output/Methylation/crossNN/crossNN_feature_coverage.tsv'
"""
model_dir = Path(args.model_dir).expanduser()
if not args.inspect:
    Path(args.predictions).parent.mkdir(parents=True, exist_ok=True)

PROBE_PREFIXES = ('cg', 'ch.', 'rs', 'ch')
#-------------------------------------------------------------------------------
# 0.3 Define functions
#-------------------------------------------------------------------------------
def as_str_list(x):
    "Coerce a label/feature container to a list of str, or return None."
    if isinstance(x, (list, tuple)) and x and all(isinstance(i, str) for i in x):
        return list(x)
    if isinstance(x, (np.ndarray, pd.Index, pd.Series)):
        arr = np.asarray(x).ravel()
        if arr.size and arr.dtype.kind in ('U', 'S', 'O'):
            return [str(i) for i in arr]
    return None


def as_matrix(x):
    "Coerce a weight container to a 2D float ndarray, or return None."
    # sklearn MLP-style estimators store a list of per-layer weight arrays
    if isinstance(x, (list, tuple)) and len(x) == 1:
        x = x[0]
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    if hasattr(x, 'detach'):                      # torch tensor, without import
        x = x.detach().cpu().numpy()
    if isinstance(x, np.ndarray) and x.ndim == 2 and x.dtype.kind == 'f':
        return x.astype(np.float32)
    return None


def as_vector(x):
    "Coerce a bias/intercept container to a 1D float ndarray, or return None."
    if isinstance(x, (list, tuple)) and len(x) == 1:
        x = x[0]
    if hasattr(x, 'detach'):
        x = x.detach().cpu().numpy()
    if isinstance(x, np.ndarray) and x.ndim == 1 and x.dtype.kind == 'f':
        return x.astype(np.float32)
    return None


def looks_like_probes(labels):
    "True if a string vector looks like Illumina probe IDs."
    if labels is None or len(labels) < 500:
        return False
    head = labels[:500]
    hits = sum(s.startswith(PROBE_PREFIXES) for s in head)
    return hits > 0.8 * len(head)


def walk(obj, path = 'root', depth = 0, max_depth = 8, seen = None):
    """Yield (path, value) for every node in a nested container/object graph.

    The crossNN pickle is a bare list holding, at unpredictable positions, the
    weight matrix, a fitted LabelEncoder and the CpG index. Names are useless
    here, so we enumerate everything and identify nodes by shape and content.
    """
    if seen is None:
        seen = set()
    if id(obj) in seen or depth > max_depth:
        return
    seen.add(id(obj))
    yield path, obj

    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f"{path}[{k!r}]", depth + 1, max_depth, seen)
    elif isinstance(obj, (list, tuple)) and not isinstance(obj, (str, bytes)):
        # Long homogeneous sequences are payloads, not containers; don't descend
        if len(obj) <= 32:
            for i, v in enumerate(obj):
                yield from walk(v, f"{path}[{i}]", depth + 1, max_depth, seen)
    elif isinstance(obj, (pd.DataFrame, pd.Series)):
        # Probe IDs live in the index or the single column, never in __dict__;
        # descending via vars() would only surface _mgr/_flags internals.
        yield f"{path}.index", obj.index
        if isinstance(obj, pd.DataFrame):
            for c in list(obj.columns)[:8]:
                yield f"{path}[{c!r}]", obj[c]
    elif hasattr(obj, '__dict__'):
        for k, v in vars(obj).items():
            yield from walk(v, f"{path}.{k}", depth + 1, max_depth, seen)


def inventory(obj):
    "Structural dump of the object graph, for --inspect and error messages."
    lines = []
    for path, v in walk(obj):
        if isinstance(v, np.ndarray):
            desc = f"ndarray shape={v.shape} dtype={v.dtype}"
        elif hasattr(v, 'detach'):
            desc = f"tensor shape={tuple(v.shape)}"
        elif isinstance(v, (list, tuple)):
            desc = f"{type(v).__name__} len={len(v)}"
        elif isinstance(v, dict):
            desc = f"dict keys={list(v)[:8]}"
        elif isinstance(v, pd.DataFrame):
            desc = f"DataFrame shape={v.shape}"
        elif isinstance(v, (pd.Index, pd.Series)):
            desc = f"{type(v).__name__} len={len(v)}"
        else:
            desc = f"{type(v).__module__}.{type(v).__name__}"
        lines.append(f"  {path}: {desc}")
    return lines


def load_artifacts(weights_path, classes_file = None, features_file = None):
    """Return (W, bias, classes, features) from the bundled crossNN pickle.

    W        : (n_classes, n_features) float32 ndarray
    bias     : (n_classes,) float32 ndarray, zeros if the model has no intercept
    classes  : methylation class labels, len == W.shape[0]
    features : CpG probe IDs, len == W.shape[1]
    """
    wpath = Path(weights_path).expanduser()
    if not wpath.exists():
        raise FileNotFoundError(f"Model pickle not found at {wpath}")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')           # sklearn version mismatch
        with open(wpath, 'rb') as fh:
            obj = pickle.load(fh)

    print(f"[crossNN] unpickled {type(obj).__module__}.{type(obj).__name__}",
          flush = True)
    if args.inspect:
        print("\n".join(inventory(obj)), flush = True)
        sys.exit(0)

    # 1) Collect candidates by shape/content across the whole graph
    mats, vecs, strs = {}, {}, {}
    for path, v in walk(obj):
        m = as_matrix(v)
        if m is not None:
            mats[path] = m
            continue
        s = as_str_list(v)
        if s is not None:
            strs[path] = s
            continue
        b = as_vector(v)
        if b is not None:
            vecs[path] = b

    feature_cands = {p: s for p, s in strs.items() if looks_like_probes(s)}
    class_cands = {p: s for p, s in strs.items()
                   if not looks_like_probes(s) and 2 <= len(s) <= 2000}

    # 2) Explicit sidecar files win over anything found in the pickle
    if features_file:
        feature_cands = {'-f': [l.strip() for l in
                                Path(features_file).read_text().splitlines()
                                if l.strip()]}
    if classes_file:
        class_cands = {'-c': [l.strip() for l in
                              Path(classes_file).read_text().splitlines()
                              if l.strip()]}

    # 3) Resolve by shape agreement: the weight matrix pins down which class and
    #    feature vectors are the real ones.
    for wpath_, W in sorted(mats.items(), key = lambda kv: -kv[1].size):
        for cpath, cls in class_cands.items():
            for fpath, feat in feature_cands.items():
                if W.shape == (len(feat), len(cls)):
                    W = W.T
                elif W.shape != (len(cls), len(feat)):
                    continue
                bias = next((b for b in vecs.values() if b.size == len(cls)),
                            np.zeros(len(cls), dtype = np.float32))
                print(f"[crossNN] weights={wpath_} classes={cpath} "
                      f"features={fpath} bias={'yes' if bias.any() else 'none'}",
                      flush = True)
                return W, bias, list(cls), list(feat)

    raise ValueError(
        f"Could not reconcile weights/classes/features in {wpath}.\n"
        f"matrices:        {[(p, m.shape) for p, m in mats.items()]}\n"
        f"probe-like:      {[(p, len(s)) for p, s in feature_cands.items()]}\n"
        f"label-like:      {[(p, len(s)) for p, s in class_cands.items()]}\n"
        f"Full structure:\n" + "\n".join(inventory(obj)))


def forward(betas, mask, W, bias, threshold):
    """Score one sample against all classes.

    Mirrors NN_classifier.predict() in the upstream NN_model.py.

    Input encoding is ternary, not beta values (training.py preprocessing):
        +1  methylated    beta >  threshold
        -1  unmethylated  beta <= threshold
         0  not observed
    The zero matches upstream's fillna(0) for probes absent from the sample, and
    matches mask_input() during training, so logits are summed without rescaling.

    The logits are then z-scored ACROSS CLASSES before the softmax:
        softmax((y - mean(y)) / std(y))
    using the population standard deviation (torch's unbiased=False). This is
    the "normalization function" referred to in the paper's Methods. Without it
    the logits are large enough that the softmax saturates at 1.0 and the
    published confidence cutoffs are meaningless.

    The model is bias-free (nn.Linear(..., bias=False)); `bias` is retained only
    so a future release carrying an intercept does not silently drop it.
    """
    n_obs = int(mask.sum())
    if n_obs == 0:
        return np.full(W.shape[0], np.nan, dtype = np.float32)
    x = np.where(mask, np.where(betas > threshold, 1.0, -1.0), 0.0).astype(np.float32)
    logits = W @ x + bias
    sd = logits.std()                      # population sd, matches unbiased=False
    if sd == 0:
        return np.full(W.shape[0], 1.0 / W.shape[0], dtype = np.float32)
    z = (logits - logits.mean()) / sd
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()
#-------------------------------------------------------------------------------
# 1.1 Load crossNN model
#-------------------------------------------------------------------------------
# The pickle may reference classes defined in the repo, so put it on sys.path
for p in (model_dir, model_dir / 'crossNN', model_dir / 'src',
          model_dir / 'scripts'):
    if p.is_dir():
        sys.path.insert(0, str(p))

W, bias, classes, features = load_artifacts(args.weights, args.classes,
                                            args.features)
print(f"[crossNN] loaded {Path(args.weights).name}: {len(classes)} classes x "
      f"{len(features)} CpGs", flush = True)
#-------------------------------------------------------------------------------
# 1.2 Read methylation data
#-------------------------------------------------------------------------------
adata = anndata.read_h5ad(args.adata)
mat = adata.layers[args.layer] if args.layer else adata.X
if hasattr(mat, 'toarray'):
    mat = mat.toarray()
mat = np.asarray(mat, dtype = np.float32)

samples = list(adata.obs_names)
probes = pd.Index(adata.var_names)

# Sanity check: crossNN expects beta values in [0,1], not M-values
finite = mat[np.isfinite(mat)]
if finite.size and (finite.min() < -0.01 or finite.max() > 1.01):
    raise ValueError(
        f"Values outside [0,1] (observed {finite.min():.2f}..{finite.max():.2f}). "
        f"crossNN expects beta values -- you are probably pointing -l at M-values.")
#-------------------------------------------------------------------------------
# 1.3 Harmonise probe IDs
#-------------------------------------------------------------------------------
# EPICv2 appends a replicate suffix (cg00000029_TC21) that will not match the
# 450k-derived model feature space. Strip it and collapse duplicates by mean.
if probes.str.contains('_').any():
    stripped = probes.str.replace(r'_.*$', '', regex = True)
    if stripped.duplicated().any():
        print(f"[crossNN] collapsing {int(stripped.duplicated().sum())} "
              f"replicate EPICv2 probes by mean", flush = True)
        df = pd.DataFrame(mat.T, index = stripped).groupby(level = 0).mean()
        mat = df.to_numpy(dtype = np.float32).T
        probes = df.index
    else:
        probes = pd.Index(stripped)
#-------------------------------------------------------------------------------
# 1.4 Project onto the model feature space
#-------------------------------------------------------------------------------
pos = probes.get_indexer(pd.Index(features))   # -1 where the probe is absent
present = pos >= 0

X = np.full((len(samples), len(features)), np.nan, dtype = np.float32)
X[:, present] = mat[:, pos[present]]

print(f"[crossNN] {present.sum()}/{len(features)} model CpGs present on this "
      f"platform", flush = True)
if present.sum() == 0:
    raise ValueError(
        "No model CpGs matched adata.var_names. Check that the h5ad is in a "
        "cg-style probe space (450k/EPIC) rather than genomic coordinates.")
#-------------------------------------------------------------------------------
# 2.1 Predict
#-------------------------------------------------------------------------------
scores = np.zeros((len(samples), len(classes)), dtype = np.float32)
n_used = np.zeros(len(samples), dtype = int)

for i, sample in enumerate(samples):
    betas = X[i]
    mask = np.isfinite(betas)          # per-sample: drops pOOBAH-masked NAs
    n_used[i] = int(mask.sum())

    if n_used[i] < args.min_features:
        warnings.warn(f"{sample}: only {n_used[i]} usable CpGs "
                      f"(< {args.min_features}); flagged as unreliable.")

    scores[i] = forward(betas, mask, W, bias, args.threshold)

score_df = pd.DataFrame(scores, index = samples, columns = classes)
#-------------------------------------------------------------------------------
# 2.2 Assemble per-sample top call
#-------------------------------------------------------------------------------
top_idx = np.nanargmax(scores, axis = 1)
pred = pd.DataFrame({
    'sample':          samples,
    'predicted_class': [classes[j] for j in top_idx],
    'score':           scores[np.arange(len(samples)), top_idx],
    'runner_up_class': [classes[j] for j in np.argsort(-scores, axis = 1)[:, 1]],
    'runner_up_score': np.sort(scores, axis = 1)[:, -2],
    'n_cpgs_used':     n_used})

pred['interpretation'] = np.where(
    pred['n_cpgs_used'] < args.min_features, 'insufficient_coverage',
    np.where(pred['score'] >= args.score_cutoff, 'match', 'no_match'))
#-------------------------------------------------------------------------------
# 3.1 Write output
#-------------------------------------------------------------------------------
pred.to_csv(args.predictions, sep = '\t', index = False)
score_df.to_csv(args.scores, sep = '\t')

pd.DataFrame({
    'sample':             samples,
    'n_model_features':   len(features),
    'n_present_platform': int(present.sum()),
    'n_cpgs_used':        n_used,
    'fraction_used':      n_used / len(features)}).to_csv(
        args.coverage, sep = '\t', index = False)

n_match = int((pred['interpretation'] == 'match').sum())
print(f"[crossNN] done: {n_match}/{len(samples)} samples above cutoff "
      f"{args.score_cutoff}", flush = True)
