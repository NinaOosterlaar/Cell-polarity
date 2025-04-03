import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def open_excel_data(file_path):
    """Function to read excel file with multiple sheets"""
    xls = pd.ExcelFile(file_path)
    sheet_names = xls.sheet_names
    data_dict = {sheet: pd.read_excel(xls, sheet) for sheet in sheet_names}
    return data_dict

def hill_equation(x, k, n):
    """Hill equation"""
    return x**n / (k**n + x**n)

def supData2(data):
    # Load data
    data.columns = ["Sample concentration", "WT_MW", "L155R_MW"]

    # Remove NaN values for fitting
    data = data.dropna(subset=["WT_MW"])

    # Define known values
    monomer_MW = 47.3  # kDa
    trimer_MW = monomer_MW * 3  # Expected trimeric MW

    # Compute fraction of trimerized protein (theta)
    data["Theta"] = data["WT_MW"] / trimer_MW

    # Fit Hill equation to data
    popt, pcov = curve_fit(hill_equation, data["Sample concentration"], data["Theta"])

    # Extract fitted k
    k_fit = popt[0]
    n_fit = popt[1]

    # Print results
    print(f"Fitted k: {k_fit:.3f}, Fitted n: {n_fit:.3f}")
    
    # Generate fitted curve
    x_vals = np.linspace(min(data["Sample concentration"]), max(data["Sample concentration"]), 100)
    theta_fit = hill_equation(x_vals, k_fit, n_fit)

    # Plot experimental data and fitted Hill equation
    plt.figure(figsize=(8, 5))
    plt.scatter(data["Sample concentration"], data["Theta"], color="blue", label="Experimental (WT)")
    plt.plot(x_vals, theta_fit, color="red", linestyle="--", label=f"Hill Fit k={k_fit:.2f}, n={n_fit:.2f}")
    plt.xlabel("Sample Concentration")
    plt.ylabel("Theta (Fraction Bound)")
    plt.title("Hill Equation Fit to Trimerization Data")
    plt.legend()
    plt.show()
    

def plot_theta_n(data):
    genotypes = ['WT', 'C49S', 'L155R']
    colors = {'WT': 'green', 'C49S': 'orange', 'L155R': 'purple'}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for genotype in genotypes:
        subset = data[data['Genotype'] == genotype]
        membrane_polarized = subset.loc[subset['State'] == 'polarized', 'Membrane concentration(dorsal) (a.u.)']
        cytosol_polarized = subset.loc[subset['State'] == 'polarized', 'Cytosol conconcentration(a.u.)']
        total_concentration_polarized = membrane_polarized + cytosol_polarized
        theta_polarized = membrane_polarized / total_concentration_polarized

        popt, pcov = curve_fit(hill_equation, cytosol_polarized, theta_polarized)
        k_fit = popt[0]
        n_fit = popt[1]

        k_error = np.sqrt(pcov[0, 0])
        n_error = np.sqrt(pcov[1, 1])

        # Compute residuals (differences between actual and fitted values)
        residuals = theta_polarized - hill_equation(cytosol_polarized, k_fit, n_fit)

        # Compute least-squares error (sum of squared residuals)
        N = len(theta_polarized)
        rmse = np.sqrt(np.sum(residuals**2)/N)

        print(f"Genotype: {genotype}, n={n_fit:.2f} ± {n_error:.2f}, k={k_fit:.2f} ± {k_error:.2f}, RMSE={rmse:.4f}")

        x_vals = cytosol_polarized.sort_values()
        theta_fit = hill_equation(x_vals, k_fit, n_fit)

        axes[0].scatter(cytosol_polarized, theta_polarized, label=f'Polarized {genotype}', color=colors[genotype])
        axes[0].plot(x_vals, theta_fit, label=f'Hill fit {genotype} (n={n_fit:.2f}, k={k_fit:.2f})', color=colors[genotype])

        axes[1].scatter(cytosol_polarized, theta_polarized, label=f'Polarized {genotype}', color=colors[genotype])
        axes[1].plot(x_vals, theta_fit, label=f'Hill fit {genotype} (n={n_fit:.2f}, k={k_fit:.2f})', color=colors[genotype])

    axes[0].set_xlabel('Cytosol concentration (a.u.)')
    axes[0].set_ylabel('Theta')
    axes[0].set_title(f'Theta for Polarized by Genotype')
    axes[0].legend()

    axes[1].set_xlabel('Cytosol concentration (a.u.) (Log scale)')
    axes[1].set_ylabel('Theta')
    axes[1].set_xscale('log')
    axes[1].set_title(f'Theta for Polarized by Genotype (Log Scale)')
    axes[1].legend()

    plt.tight_layout()
    plt.show()

def plot_theta_fixed_n(data, n):
    genotypes = ['WT', 'C49S', 'L155R']
    colors = {'WT': 'green', 'C49S': 'orange', 'L155R': 'purple'}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    for genotype in genotypes:
        subset = data[data['Genotype'] == genotype]

        membrane_polarized = subset.loc[subset['State'] == 'polarized', 'Membrane concentration(dorsal) (a.u.)']
        cytosol_polarized = subset.loc[subset['State'] == 'polarized', 'Cytosol conconcentration(a.u.)']
        total_concentration_polarized = membrane_polarized + cytosol_polarized
        theta_polarized = membrane_polarized / total_concentration_polarized

        # Fit only k while keeping n fixed
        popt, pcov = curve_fit(lambda x, k: hill_equation(x, k, n), cytosol_polarized, theta_polarized, p0=[1.0])
        k_fit = popt[0]
        k_error = np.sqrt(pcov[0, 0])

        # Compute residuals (differences between actual and fitted values)
        residuals = theta_polarized - hill_equation(cytosol_polarized, k_fit, n)

        # Compute least-squares error (sum of squared residuals)
        N = len(theta_polarized)
        rmse = np.sqrt(np.sum(residuals**2)/N)

        print(f"Genotype: {genotype}, n={n}, k={k_fit:.2f} ± {k_error:.2f}, RMSE={rmse:.4f}")

        # Sort x values for smooth plotting
        x_vals = np.linspace(cytosol_polarized.min(), cytosol_polarized.max(), 100)
        theta_fit = hill_equation(x_vals, k_fit, n)

        # Scatter plot + fitted curve
        axes[0].scatter(cytosol_polarized, theta_polarized, label=f'Polarized {genotype}', color=colors[genotype])
        axes[0].plot(x_vals, theta_fit, label=f'Hill fit {genotype} (n={n}, k={k_fit:.2f})', color=colors[genotype])

        axes[1].scatter(cytosol_polarized, theta_polarized, label=f'Polarized {genotype}', color=colors[genotype])
        axes[1].plot(x_vals, theta_fit, label=f'Hill fit {genotype} (n={n}, k={k_fit:.2f})', color=colors[genotype])

    # Set plot labels and legends
    axes[0].set_xlabel('Cytosol concentration (a.u.)')
    axes[0].set_ylabel('Theta')
    axes[0].set_title(f'Theta for Polarized by Genotype')
    axes[0].legend()

    axes[1].set_xlabel('Cytosol concentration (a.u.) (Log scale)')
    axes[1].set_ylabel('Theta')
    axes[1].set_xscale('log')
    axes[1].set_title(f'Theta for Polarized by Genotype (Log Scale)')
    axes[1].legend()

    plt.tight_layout()
    plt.show()


# Load data
data = pd.read_excel("DataSet_CellPolarity.xlsx", sheet_name="MAIN data")
plot_theta_n(data)
data = pd.read_excel("DataSet_CellPolarity.xlsx", sheet_name="MAIN data")
plot_theta_fixed_n(data, n=0.5)
plot_theta_fixed_n(data, n=1.0)
plot_theta_fixed_n(data, n=1.5)
plot_theta_fixed_n(data, n=2.0)
data = open_excel_data("DataSet_CellPolarity.xlsx")["SupData2"]
supData2(data)

