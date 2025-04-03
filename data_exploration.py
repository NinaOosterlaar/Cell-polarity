import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def open_excel_data(file_path):
    """Function to read excel file with multiple sheets"""
    # Read the excel file
    xls = pd.ExcelFile(file_path)
    # Get the sheet names
    sheet_names = xls.sheet_names
    # Create a dictionary to store the data
    data_dict = {}
    # Loop through the sheet names
    for sheet in sheet_names:
        # Read the data
        data = pd.read_excel(xls, sheet)
        # Store the data in the dictionary
        data_dict[sheet] = data
    return data_dict


def visualize_boxplots_states(data):
    """Creates four different plots for each genotype, where each plot has 4 subplots 
        showing the concentration for polarized and non-polarized states"""
    df = data['MAIN data']
    # Extract columns that contain 'concentration' in their name
    concentration_columns = [col for col in df.columns if 'concentration' in col.lower()]
    genotypes = df['Genotype'].unique()

    for genotype in genotypes:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))  # 2x2 subplots
        df_genotype = df[df['Genotype'] == genotype]

        for i, column in enumerate(concentration_columns):
            row, col = divmod(i, 2)
            sns.boxplot(
                data=df_genotype[df_genotype["State"].isin(["polarized", "non-polarized"])],
                x="State",
                y=column,
                ax=axes[row, col]
            )
            axes[row, col].set_title(f"{column} - {genotype}")
            axes[row, col].set_xlabel("State")
            axes[row, col].set_ylabel("Concentration")

        plt.tight_layout()
        plt.show()
        
        
def visualize_boxplots_genotypes(data):
    """Creates four different plots for each concentration, where each plot has two subplots 
        showing the concentration for each genotype for polarized and non-polarized states"""
    df = data['MAIN data']
    # Extract columns that contain 'concentration' in their name
    concentration_columns = [col for col in df.columns if 'concentration' in col.lower()]
    genotypes_order = ["WT", "L155R", "C49S"]  # Define the order of genotypes
    
    for concentration in concentration_columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # 1x2 subplots
        
        for i, state in enumerate(["polarized", "non-polarized"]):
            sns.boxplot(
                data=df[df["State"] == state],
                x="Genotype",
                y=concentration,
                ax=axes[i],
                palette="Set2",  # Use a different color palette
                order=genotypes_order  # Ensure the order of genotypes
            )
            axes[i].set_title(f"{concentration} - {state}", fontsize=14)
            axes[i].set_xlabel("Genotype", fontsize=12)
            axes[i].set_ylabel("Concentration", fontsize=12)
            axes[i].tick_params(axis='both', which='major', labelsize=10)
        
        plt.tight_layout()
        plt.show()
        
        
def heatmap(data):
    """Create two heatmaps plotting the concentration values for each genotype and state."""
    df = data['MAIN data']
    
    # Extract concentration columns (assumes they contain 'concentration' in their name)
    concentration_columns = [col for col in df.columns if 'concentration' in col.lower()]
    
    # Create one heatmap for polarized and one for non-polarized state
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # 1x2 layout

    for i, state in enumerate(["polarized", "non-polarized"]):
        subset = df[df["State"] == state]
        
        # Create a pivot table: Rows = Genotype, Columns = Concentration types, Values = Concentration
        heatmap_data = subset.pivot_table(index="Genotype", values=concentration_columns, aggfunc="mean")
        
        # Rename columns for shorter names
        heatmap_data.columns = [
            'cytosol' if 'Cytosol' in col else 
            'total' if 'total' in col else 
            'dorsal' if 'dorsal' in col else 
            'ventral' if 'ventral' in col else col
            for col in heatmap_data.columns
        ]

        sns.heatmap(
            data=heatmap_data, 
            ax=axes[i], 
            cmap="viridis", 
            annot=True, 
            fmt=".2f"
        )
        axes[i].set_title(f"{state.capitalize()} State")
        axes[i].set_xlabel("Concentration Type")
        axes[i].set_ylabel("Genotype")

    plt.tight_layout()
    plt.show()
            

def scatter_protein_size(data):
    df = data['SupData2']
    sample_concentration = df['Sample concentration']
    wildtype_protein_weight = df['WT_Average Molecular weight (kDa)']
    L155R_protein_weight = df['L155R_Average Molecular Weight(kDa)']
    
    plt.plot(sample_concentration, wildtype_protein_weight, 'o', label='Wildtype')
    plt.plot(sample_concentration, L155R_protein_weight, 'o', label='L155R')
    plt.xlabel('Sample concentration')
    plt.ylabel('Average Molecular weight (kDa)')
    plt.title('Protein size vs Sample concentration')
    plt.legend()
    plt.show()
    
    
def another_boxplot(data):
    df1 = data['MAIN data']
    df2 = data['SupData3']
    
    dictionary = {"WT": {"Cytosol": [], "Total": [], "Dorsal": [], "Ventral": []},
                  "L155R": {"Cytosol": [], "Total": [], "Dorsal": [], "Ventral": []},
                  "C49S": {"Cytosol": [], "Total": [], "Dorsal": [], "Ventral": []},
                  "Other Protein": {"Cytosol": [], "Total": [], "Dorsal": [], "Ventral": []}}
    
    df1_wt = df1[df1['Genotype'] == 'WT']
    df1_l155r = df1[df1['Genotype'] == 'L155R']
    df1_c49s = df1[df1['Genotype'] == 'C49S']
    
    for column in df1.columns:
        if "Cytosol" in column:
            dictionary["WT"]["Cytosol"].extend(df1_wt[column].values)
            dictionary["L155R"]["Cytosol"].extend(df1_l155r[column].values)
            dictionary["C49S"]["Cytosol"].extend(df1_c49s[column].values)
        elif "total" in column:
            dictionary["WT"]["Total"].extend(df1_wt[column].values)
            dictionary["L155R"]["Total"].extend(df1_l155r[column].values)
            dictionary["C49S"]["Total"].extend(df1_c49s[column].values)
        elif "dorsal" in column:
            dictionary["WT"]["Dorsal"].extend(df1_wt[column].values)
            dictionary["L155R"]["Dorsal"].extend(df1_l155r[column].values)
            dictionary["C49S"]["Dorsal"].extend(df1_c49s[column].values)
        elif "ventral" in column:
            dictionary["WT"]["Ventral"].extend(df1_wt[column].values)
            dictionary["L155R"]["Ventral"].extend(df1_l155r[column].values)
            dictionary["C49S"]["Ventral"].extend(df1_c49s[column].values)
    
    for column in df2.columns:
        if "Cytosol" in column:
            dictionary["Other Protein"]["Cytosol"].extend(df2[column].values)
        elif "total" in column:
            dictionary["Other Protein"]["Total"].extend(df2[column].values)
        elif "dorsal" in column:
            dictionary["Other Protein"]["Dorsal"].extend(df2[column].values)
        elif "ventral" in column:
            dictionary["Other Protein"]["Ventral"].extend(df2[column].values)
            
    print(dictionary["Other Protein"]["Dorsal"])
            
    # Ensure equal length for each list in the category for boxplot creation
    categories = ["Cytosol", "Total", "Dorsal", "Ventral"]
    
    # Create a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))  # Adjust size as needed
    axes = axes.flatten()  # Flatten to access each subplot easily
    
    for idx, category in enumerate(categories):
        max_len = max(len(dictionary[genotype][category]) for genotype in dictionary)
        
        for genotype in dictionary:
            # Extend each list to the maximum length by adding NaN values for missing data
            while len(dictionary[genotype][category]) < max_len:
                dictionary[genotype][category].append(np.nan)
                
        data_to_plot = {
            "WT": dictionary["WT"][category],
            "L155R": dictionary["L155R"][category],
            "C49S": dictionary["C49S"][category],
            "Other Protein": dictionary["Other Protein"][category]
        }
        
        # Convert to DataFrame for seaborn compatibility
        plot_data = pd.DataFrame(data_to_plot)
        
        # Create the boxplot in the corresponding subplot
        sns.boxplot(data=plot_data, ax=axes[idx])
        
        # Set plot labels and title
        axes[idx].set_title(f"Boxplot of {category}")
        axes[idx].set_xlabel("Genotype")
        axes[idx].set_ylabel("Values")
    
    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.show()
    

def average_concentrations(data):
    data = data['MAIN data']
    genotypes = data['Genotype'].unique()
    states = data['State'].unique()
    for genotype in genotypes:
        for state in states:
            subset = data[(data['Genotype'] == genotype) & (data['State'] == state)]
            dorsal_concentration = subset['Membrane concentration(dorsal) (a.u.)']
            cytosol_concentration = subset['Cytosol conconcentration(a.u.)']
            average_dorsal_concentration = dorsal_concentration.mean()
            average_cytosol_concentration = cytosol_concentration.mean()
            print(f"Genotype: {genotype}, State: {state}")
            print(f"Average Dorsal Concentration: {average_dorsal_concentration}")
            print(f"Average Cytosol Concentration: {average_cytosol_concentration}")
            print()
    
    
    
    
if __name__ == "__main__":
    data = open_excel_data("DataSet_CellPolarity.xlsx")
    # visualize_boxplots_states(data)
    # visualize_boxplots_genotypes(data)
    # heatmap(data)
    # scatter_protein_size(data)
    # another_boxplot(data)
    average_concentrations(data)
    