import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from scipy import stats
import time
from scipy.stats import binned_statistic_2d
import os

os.makedirs("figures", exist_ok=True)

# Mock genotype matrix, assuming first 4034 SNPs are on chromosome a, rest is chromsome b. Needs to be replaced with real genotypes.
df = np.random.randint(0, 3, size=(10000, 8170))

# Mock covariates. Needs to be replaced with real covariates.
covs = np.random.randn(10000, 7)

# Mean center SNPs
df1 = np.zeros_like(df)
df1 = df1.astype(float)
for i in range(df.shape[1]):
    df1[:,i] = df[:,i] - np.mean(df[:,i])

# Function to sample two-way interacting SNPs
def sample_column_pairs(n_cols, n_pairs):
    pairs = set()
    while len(pairs) < n_pairs:
        # pick two distinct columns
        i, j = np.random.choice(n_cols, size=2, replace=False)
        # store as an ordered tuple so (i, j) == (j, i)
        pairs.add(tuple(sorted((i, j))))
    return np.array(list(pairs))


# Function to obtain SNP indices exlucing SNPs from a set
def sample_excluding(low, high, exclude, size=1):
    """Sample random integers in [low, high) excluding given numbers."""
    exclude = set(exclude)
    candidates = [x for x in range(low, high) if x not in exclude]
    temp = np.random.choice(candidates, size=size, replace=False)
    return [temp]

n, m = df1.shape

# Build covariate design with intercept
C = np.column_stack([np.ones(n), np.array(covs)])   # shape (n, p), here p = 8

# Precompute projection onto covariates
CtC_inv = np.linalg.inv(C.T @ C)
B = CtC_inv @ C.T

# Residualize all SNPs at once
df1_res = df1 - C @ (B @ df1)

# Lambda proportions of the variance explained by the interaction term
lambd_scales = np.linspace(0.001,0.171,100)

chromosome_a_size = 4034

# Lists to hold parameters
rho_list = []
lambd_list = []
p_val = []
t_stat = []


start = time.time()


# Loop for simulating phenotypes and performing OLS. Note that OLS is performed here as a quick
for seed in range(len(lambd_scales)):
    print(seed)
    lambd = lambd_scales[seed]

    random.seed(seed)
    np.random.seed(seed)  

    # Get SNP pairs from chromosome a
    pairs = sample_column_pairs(chromosome_a_size, 100)

    # Vectorized Z construction
    A1 = df1[:, pairs[:, 0]]
    B1 = df1[:, pairs[:, 1]]
    Z = A1 * B1
    Z -= Z.mean(axis=0, keepdims=True) # center each product column

    # SNP indices for additve SNPs on chromosome b
    G_ix = sample_excluding(chromosome_a_size, df1.shape[1], exclude=pairs.flatten(), size=1000)[0]

    # Coefficients for additive SNPs
    G_betas = np.random.normal(loc=0, scale=1, size=len(G_ix))

    # Coefficients for interaction features
    Z_betas = np.random.normal(loc=0, scale=1, size=Z.shape[1])
    
    G_realized = (G_betas * df1[:,G_ix])
    Z_realized = (Z_betas * Z)

    # Realized vector of additive SNP effects
    epsilon_kappa = np.sum(G_realized, axis=1)

    # Realized vector of interaction effects
    u = np.sum(Z_realized, axis=1)    

    # Realized vector of noise
    epsilon_e = np.random.normal(loc=0, scale=1, size=len(G_realized))

    # 
    epsilon_kappa = epsilon_kappa/np.sqrt(np.var(epsilon_kappa) + np.var(epsilon_e))
    u = u/np.sqrt(np.var(epsilon_kappa) + np.var(epsilon_e))
    epsilon_e = epsilon_e/np.sqrt(np.var(epsilon_kappa) + np.var(epsilon_e))

    # Var(u)
    var_u = lambd/(1-lambd)

    # 
    target_std = np.sqrt(var_u)
    
    std_x = np.std(u, ddof=1)
    u = u * (target_std / std_x)

    # Phenotype y
    y = u + epsilon_kappa + epsilon_e

    # Sample target SNP indices on chromosome a excluding SNPs in the interaction term
    noncausal_ix = sample_excluding(0, chromosome_a_size, exclude=np.concatenate([pairs.flatten(), G_ix]), size=3834)[0]

    X = df1[:, noncausal_ix]

    # Residualize to remove covariates
    X_res = df1_res[:, noncausal_ix]
    y_res  = y  - C @ (B @ y)

    # Precompute sufficient stats
    xTx   = np.sum(X_res * X_res, axis=0)
    XTy   = X_res.T @ y_res
    ySS   = y_res  @ y_res
    
    # degrees of freedom: n - (#covariates incl. intercept) - 1 (the SNP)
    p = C.shape[1]
    df_resid = n - p - 1

    # Obtain coefficients
    beta  = XTy  / xTx
    SSE  = ySS  - (XTy**2)  / xTx
    MSE  = SSE  / df_resid

    # Obtain SE
    se  = np.sqrt(MSE  / xTx)

    # Obtain t-statistics
    tvals  = beta  / se
    
    # Obtain p-values
    pvals  = 2.0 * stats.t.sf(np.abs(tvals),  df_resid)
    
    # Add p-values and t-statistics to the lists
    p_val.extend(pvals.tolist())
    t_stat.extend(tvals.tolist())
    
    # rho correlations corr(u, x_j) for all j
    u_m  = u - u.mean()
    X_m  = X - X.mean(axis=0, keepdims=True)
    corr_vec = (u_m @ X_m) / ((n - 1) * u_m.std(ddof=1) * X_m.std(axis=0, ddof=1))

    # Add rho and lambda values to the list
    rho_list.extend(corr_vec.tolist())
    lambd_list.extend([lambd] * X.shape[1])
    

rho_list = np.array(rho_list)
lambd_list = np.array(lambd_list)
p_val = np.array(p_val)
t_stat = np.array(t_stat)


end = time.time()

print(f"time elapsed: {end-start}")

x = rho_list
y = lambd_list
z = np.abs(t_stat)

# Choose grid resolution 
nx = ny = 50
level = 5.45  # GWAS significance threshold for highlighting
n = df.shape[0]


# Aggregate z per (x,y) bin — mean or median
stat, xedges, yedges, _ = binned_statistic_2d(
    x, y, z, statistic='mean', bins=[nx, ny])

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
pcm = ax.pcolormesh(xedges, yedges, stat.T, cmap='inferno', shading='auto')
fig.colorbar(pcm, ax=ax, label=r'abs($z$)')

# Create boolean masks for cells above the threshold
mask = np.isfinite(stat) & (stat >= level)

# Draw borders for a mask with a specific color
def draw_borders(mask, color, linewidth=1.2):
    for i in range(nx):
        for j in range(ny):
            if mask[i, j]:
                # determine which edges to draw
                draw_left = (i == 0 or not mask[i - 1, j])
                draw_right = (i == nx - 1 or not mask[i + 1, j])
                draw_bottom = (j == 0 or not mask[i, j - 1])
                draw_top = (j == ny - 1 or not mask[i, j + 1])
                
                # draw individual edge segments
                x0, x1 = xedges[i], xedges[i + 1]
                y0, y1 = yedges[j], yedges[j + 1]
                
                if draw_left:
                    ax.plot([x0, x0], [y0, y1], color=color, linewidth=linewidth)
                if draw_right:
                    ax.plot([x1, x1], [y0, y1], color=color, linewidth=linewidth)
                if draw_bottom:
                    ax.plot([x0, x1], [y0, y0], color=color, linewidth=linewidth)
                if draw_top:
                    ax.plot([x0, x1], [y1, y1], color=color, linewidth=linewidth)

# Draw borders for high values (cyan)
draw_borders(mask, 'cyan')

ax.set_xlim(-0.3, 0.3)

plt.title(fr'$n={n}$')
ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$\lambda$')
plt.show()
fig.savefig(f'figures/sim_lambda_rho_n={n}_max.png', dpi=300)
