#!/usr/bin/env python3
import numpy as np
import time
import sys

for k in (40000, 60000, 80000, 100000):
    print(f"=== k={k} ===", flush=True)
    t_alloc = time.perf_counter()
    A = np.random.default_rng(0).normal(size=(k, k)).astype(np.float64)
    P = A @ A.T
    del A
    P += float(k) * np.eye(k, dtype=np.float64)
    print(f"k={k} build done in {time.perf_counter()-t_alloc:.1f}s, starting cholesky...", flush=True)
    t0 = time.perf_counter()
    try:
        L = np.linalg.cholesky(P)
        print(f"k={k} cholesky OK in {time.perf_counter()-t0:.1f}s", flush=True)
    except Exception as e:
        print(f"k={k} cholesky FAILED: {e!r}", flush=True)
        sys.exit(1)
    del P, L
print("ALL_SIZES_OK", flush=True)
