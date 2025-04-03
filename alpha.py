import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Define the cooperative binding model
def cooperative_binding(c, beta, alpha):
    return beta * c ** alpha

# Function to plot cooperativity for different genotypes
def plot_cooperativity(data):
    genotypes = ['WT', 'C49S', 'L155R']
    colors = {'WT': 'green', 'C49S': 'orange', 'L155R': 'purple'}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for genotype in genotypes:
        subset = data[data['Genotype'] == genotype]

        # Extract the polarized membrane and cytosol concentrations
        membrane_polarized = subset.loc[subset['State'] == 'polarized', 'Membrane concentration(dorsal) (a.u.)']
        cytosol_polarized = subset.loc[subset['State'] == 'polarized', 'Cytosol conconcentration(a.u.)']

        # Fit the cooperative model
        popt, pcov = curve_fit(cooperative_binding, cytosol_polarized, membrane_polarized, p0=[1.0, 1.0])
        beta_fit, alpha_fit = popt
        alpha_error = np.sqrt(np.diag(pcov))[1]  # Error for alpha

        # Compute residuals (differences between actual and fitted values)
        residuals = membrane_polarized - cooperative_binding(cytosol_polarized, *popt)

        # Compute least-squares error (sum of squared residuals)
        N = len(membrane_polarized)
        rmse = np.sqrt(np.sum(residuals**2) / N)

        print(f"Genotype: {genotype}, beta={beta_fit:.2f}, alpha={alpha_fit:.2f} ± {alpha_error:.2f}, RMSE={rmse:.4f}")

        # Sort x values for smooth plotting
        x_vals = np.linspace(cytosol_polarized.min(), cytosol_polarized.max(), 100)
        y_fit = cooperative_binding(x_vals, *popt)

        # Scatter plot + fitted curve
        axes[0].scatter(cytosol_polarized, membrane_polarized, label=f'{genotype}', color=colors[genotype])
        axes[0].plot(x_vals, y_fit, label=f'Cooperative fit {genotype} (alpha={alpha_fit:.2f})', color=colors[genotype])

        axes[1].scatter(cytosol_polarized, membrane_polarized, label=f'{genotype}', color=colors[genotype])
        axes[1].plot(x_vals, y_fit, label=f'Cooperative fit {genotype} (alpha={alpha_fit:.2f})', color=colors[genotype])

    # Set plot labels and legends
    axes[0].set_xlabel('Cytosol concentration (a.u.)')
    axes[0].set_ylabel('Membrane concentration (a.u.)')
    axes[0].set_title(f'Membrane vs Cytosol Concentration by Genotype')
    axes[0].legend()

    axes[1].set_xlabel('Cytosol concentration (a.u.) (Log scale)')
    axes[1].set_ylabel('Membrane concentration (a.u.)')
    axes[1].set_xscale('log')
    axes[1].set_title(f'Membrane vs Cytosol Concentration by Genotype (Log Scale)')
    axes[1].legend()

    plt.tight_layout()
    plt.show()

# Load data and run the function
data = pd.read_excel("DataSet_CellPolarity.xlsx", sheet_name="MAIN data")
plot_cooperativity(data)
