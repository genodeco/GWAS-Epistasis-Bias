import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import os

os.makedirs("figures", exist_ok=True)

# Mock genotype matrix, assuming first 4034 SNPs are on chromosome a, rest is chromsome b. Needs to be replaced with real genotypes.
df = np.random.randint(0, 3, size=(10000, 8170))

# Mean center SNPs
df1 = np.zeros_like(df)
df1 = df1.astype(float)
for i in range(df.shape[1]):
    df1[:,i] = df[:,i] - np.mean(df[:,i])

# Mock covariates. Needs to be replaced with real covariates.
covs = np.random.randn(10000, 7)

def sample_column_pairs(n_cols, n_pairs):
    pairs = set()
    while len(pairs) < n_pairs:
        i, j = np.random.choice(n_cols, size=2, replace=False)
        if i > j:
            i, j = j, i
        pairs.add((i, j))
    return np.array(list(pairs), dtype=int)


def sample_excluding(low, high, exclude, size=1):

    if low >= high:
        raise ValueError("low must be less than high.")

    exclude = np.asarray(exclude, dtype=int)
    # Clip exclude to valid range
    exclude = exclude[(exclude >= low) & (exclude < high)]

    mask = np.ones(high - low, dtype=bool)
    mask[exclude - low] = False
    candidates = np.nonzero(mask)[0] + low

    if size > candidates.size:
        raise ValueError(f"Requested {size} samples but only {candidates.size} available after exclusion.")
    
    return np.random.choice(candidates, size=size, replace=False)

n, m = df1.shape

# Build covariate design with intercept
C = np.column_stack([np.ones(n), np.asarray(covs)])

CtC = C.T @ C

# Residualize all SNPs at once
df1 = df1 - C @ np.linalg.solve(CtC, C.T @ df1)


n_cols = 4034 # interacting SNPs on chromosome a
n_iters = 100 # 100 different Z interaction matrices
n_g = 4034 - 200 # Same chr
#n_g = 4135 # Diff chr

batch_size = 256 # adjust for speed/memory
n_pairs_list = [100] # 100 2-way interaction features
rho_list_all = []

for n_pairs in n_pairs_list:
    print(f"n_pairs = {n_pairs}")
    rho_per_iter = []

    for _ in range(n_iters):
        print(_)

        random.seed(_)
        np.random.seed(_)
        
        # Build Z from pairwise products, mean-center each column
        pairs = sample_column_pairs(n_cols, n_pairs)
        Z = df1[:, pairs[:, 0]] * df1[:, pairs[:, 1]]
        Z -= Z.mean(axis=0, keepdims=True)  # center the product columns

        Z = Z - C @ np.linalg.solve(CtC, C.T @ Z) 

        try:
            R = np.linalg.qr(Z, mode='r')  # returns only R (NumPy ≥ 1.10)
        except TypeError:
            # Fallback if 'mode="r"' not available
            R = np.linalg.qr(Z, mode='reduced')[1]

        # Choose distinct g columns excluding anything used in Z
        exclude_cols = np.unique(pairs.flatten())

        g_indices = sample_excluding(0, 4034, exclude_cols, size=n_g) # same chr
        #g_indices = sample_excluding(4035, 8170, exclude_cols, size=n_g) # dif chr

        # Compute rhos in batches: rho = ||Q^T G|| / ||G||
        rhos_this_iter = np.empty(n_g, dtype=Z.dtype)
        start = 0
        while start < n_g:
            end = min(start + batch_size, n_g)
            G = df1[:, g_indices[start:end]]                  # (n x b)
            G_norms = np.linalg.norm(G, axis=0)               # (b,)

            # Y = Z^T G (k x b)
            Y = Z.T @ G

            # Solve R^T T = Y  -> T = R^{-T} Y  (k x b)
            T = np.linalg.solve(R.T, Y)

            # ||Q^T G|| = sqrt(sum(T^2, axis=0))
            proj_norms = np.sqrt(np.sum(T * T, axis=0))

            # safe division
            with np.errstate(divide='ignore', invalid='ignore'):
                rhos_batch = proj_norms / G_norms
                rhos_batch[~np.isfinite(rhos_batch)] = 0.0  # handles any zero-norm edge cases

            rhos_this_iter[start:end] = rhos_batch
            start = end

        rho_per_iter.append(rhos_this_iter)

        del Z

    rho_list_all.append(np.vstack(rho_per_iter))

rho_list_all = np.stack(rho_list_all, axis=0)

fig = plt.figure(figsize=(6, 6))
plt.hist(rho_list_all.flatten(), bins=100, log=True)
plt.xlabel(r'$|\rho|$')
plt.ylabel('counts (log scale)')
plt.title(fr'Distribution of $\rho_{{max}}$')

plt.axvline(x=np.mean(rho_list_all.flatten()), color='white', linestyle='--', linewidth=1)

fig.savefig('figures/rho_max_distribution.png', dpi=300)
plt.show()
