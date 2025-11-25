import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import os

os.makedirs("figures", exist_ok=True)

n = 100_000
x = 5.45

rho_vals = np.linspace(-0.3, 0.3, 400)
lambda_vals = np.linspace(0, 0.171, 400)

RHO, LAMBDA = np.meshgrid(rho_vals, lambda_vals)

# Model formulas
MU = (RHO * np.sqrt(LAMBDA * n)) / np.sqrt(1 - LAMBDA * RHO**2)
SIGMA = np.sqrt((1 - LAMBDA * RHO**2) / (1 - LAMBDA))
alpha_t = norm.cdf(-(SIGMA*(x-MU))) + norm.cdf(-(SIGMA*(x+MU)))
alpha_w = 2 * norm.cdf(-x)
R = alpha_t/alpha_w

fig = plt.figure(figsize=(8, 6))
plt.pcolormesh(RHO, LAMBDA, R, shading='auto', cmap='inferno')
plt.colorbar(label=r'$R$  (for $x=\pm5.45$)')


contour_level = 1
cs = plt.contour(RHO, LAMBDA, R,
                 levels=[contour_level], colors='cyan', linewidths=1.2)

contour_level = 0.5/alpha_w
cs = plt.contour(RHO, LAMBDA, R,
                 levels=[contour_level], colors='cyan', linewidths=1.2, linestyles='dashed')
contour_level = 0.9/alpha_w
cs = plt.contour(RHO, LAMBDA, R,
                 levels=[contour_level], colors='cyan', linewidths=1.2, linestyles='dashed')

contour_level = 0.1/alpha_w
cs = plt.contour(RHO, LAMBDA, R,
                 levels=[contour_level], colors='cyan', linewidths=1.2, linestyles='dashed')

plt.xlabel(r'$\rho$')
plt.ylabel(r'$\lambda$')
plt.title(fr'$n={n}$')
plt.show()

fig.savefig(f'figures/lambda_rho_n={n}.png', dpi=300)




lambda_val_list = [0.055, 0.119, 0.171]

for lambda_val in lambda_val_list:

    x = 5.45
    n_vals = np.linspace(10_000, 1_000_000, 400)
    rho_val_lim = 0.032
    rho_vals = np.linspace(-rho_val_lim, rho_val_lim, 400)  
    
    # Create grid
    N, RHO = np.meshgrid(n_vals, rho_vals)
    
    # Model formulas
    MU = (RHO * np.sqrt(N * lambda_val)) / np.sqrt(1 - lambda_val * RHO**2)
    SIGMA = np.sqrt((1 - lambda_val * RHO**2) / (1 - lambda_val))
    alpha_t = norm.cdf(-(SIGMA*(x-MU))) + norm.cdf(-(SIGMA*(x+MU)))
    alpha_w = 2 * norm.cdf(-x)
    R = alpha_t/alpha_w
    
    # Plot heatmap
    fig = plt.figure(figsize=(8,6))
    plt.pcolormesh(N, RHO, R, shading='auto', cmap='inferno')
    plt.colorbar(label=r'$R$  (for $x=\pm5.45$)')
    
    contour_level = 1
    cs = plt.contour(N, RHO, R, levels=[contour_level],
                     colors='cyan', linewidths=2)
    contour_level = 0.5/alpha_w
    cs = plt.contour(N, RHO, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    contour_level = 0.1/alpha_w
    cs = plt.contour(N, RHO, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    contour_level = 0.9/alpha_w
    cs = plt.contour(N, RHO, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    plt.axvline(x=210145, color='white', linestyle='--', linewidth=1)
    
    
    plt.ylabel(r'$\rho$')
    plt.xlabel(r'$n$')
    plt.title(fr'$\lambda$ = {np.round(lambda_val, 3)}')
    plt.show()
    
    fig.savefig(f'figures/rho_n_lambda={np.round(lambda_val, 3)}_rholim={np.round(rho_val_lim, 3)}.png', dpi=300)



rho_max_list = np.array([0.021 , 0.032, 0.042, 0.024, 0.035, 0.849])

for rho in rho_max_list:
    x = 5.45
    n_vals = np.linspace(10_000, 1_000_000, 200)
    lambda_vals = np.linspace(0, 0.171, 200)
    
    # Create grid
    N, LAMBDA = np.meshgrid(n_vals, lambda_vals)
    
    # Model formulas
    MU = (rho * np.sqrt(N * LAMBDA)) / np.sqrt(1 - LAMBDA * rho**2)
    SIGMA = np.sqrt((1 - LAMBDA * rho**2) / (1 - LAMBDA))
    alpha_t = norm.cdf(-(SIGMA*(x-MU))) + norm.cdf(-(SIGMA*(x+MU)))
    alpha_w = 2 * norm.cdf(-x)
    R = alpha_t/alpha_w
    
    # Plot
    fig = plt.figure(figsize=(8,6))
    plt.pcolormesh(N, LAMBDA, R, shading='auto', cmap='inferno')
    plt.colorbar(label=r'$R$  (for $x=\pm5.45$)')
    
    contour_level = 0.5/alpha_w
    cs = plt.contour(N, LAMBDA, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    contour_level = 0.9/alpha_w
    cs = plt.contour(N, LAMBDA, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    contour_level = 0.1/alpha_w
    cs = plt.contour(N, LAMBDA, R, levels=[contour_level],
                     colors='cyan', linewidths=2, linestyles='dashed')
    
    contour_level = 1.001
    cs = plt.contour(N, LAMBDA, R, levels=[contour_level],
                     colors='cyan', linewidths=2)
    
    
    plt.axvline(x=210145, color='white', linestyle='--', linewidth=1)
    

    plt.xlabel(r'$n$')
    plt.ylabel(r'$\lambda$')
    plt.title(fr'$\rho$ = {np.round(rho, 3)}')
    plt.show()
    
    fig.savefig(f'figures/lambda_n_rho={np.round(rho, 3)}.png', dpi=300)
