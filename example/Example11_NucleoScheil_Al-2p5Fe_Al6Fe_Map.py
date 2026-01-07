
import pandas as pd, numpy as np
import polars as pl
import equilipy as eq
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os

axfont= {
    'size': 16,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
} 
def process_data(df: pd.DataFrame, phase_map: dict) -> pd.DataFrame:
    """
    Loads data from a CSV and creates a 'Configuration' column based on active phases.

    Args:
        filename: The path to the CSV file.
        phase_map: A dictionary mapping phase names to short labels.

    Returns:
        A pandas DataFrame with the new 'Configuration' column.
    """
    
    # Use a small tolerance to handle floating-point inaccuracies
    tolerance = 1e-9
    phase_cols = list(phase_map.keys())
    
    # Identify which phases are active (present) in each row
    active_phases = df[phase_cols] > tolerance

    def get_config_label(row_of_bools):
        """Generates a sorted, combined label for active phases in a row."""
        # Select the phase names where the row value is True
        active_labels = [phase_map[col] for col, is_active in row_of_bools.items() if is_active]
        
        if not active_labels:
            return "None"
            
        # Sort to ensure consistency (e.g., "B+e" is the same as "e+B")
        active_labels.sort()
        return "+".join(active_labels)

    # Apply the function to each row to generate the configuration string
    df['Configuration'] = active_phases.apply(get_config_label, axis=1)
    
    return df

def plot_phase_map(df: pd.DataFrame, description_map: dict):
    """
    Creates and formats a scatter plot of the phase configurations.

    Args:
        df: The processed DataFrame containing the data to plot.
        description_map: A map from short labels back to full phase names for the legend.
    """
    # 1. Set up the plot
    fig, ax = plt.subplots(figsize=(5, 4.5))
    deep = [
    "#4C72B0",  # Blue
    "#CCB974",  # Yellow-Green
    "#55A868",  # Green
    "#C44E52",  # Red
    "#8172B3",  # Purple
    "#DD8452",  # Orange
    "#DA8BC3",  # Pink
    "#8C8C8C",  # Gray
    "#937860",  # Brown
    "#64B5CD",   # Light Blue
    "#ffa500",  # SteelBlue
    "#4169e1"
]
    
    # 2. Create the scatter plot
    # Using square markers with no edge color helps simulate a filled region map
    # tab20_dark = mpl.cm.get_cmap('tab20', 20)
    # dark_tab_palette = [tab20_dark(i) for i in range(1, 20, 2)][:12]
    sns.scatterplot(
        data=df,
        x='dT_AL13FE4',
        y='dT_Al6Fe_AL6MN(s)',
        hue='Configuration',
        palette=deep,  # A color palette with good visual distinction
        s=50,             # Marker size
        marker='s',       # Square markers fill space well
        edgecolor='none',
        ax=ax,
        legend=False
    )

    # 3. Format the plot and axes
    ax.set_xlabel(r'$\boldsymbol{\Delta}$ T$^{\boldsymbol{*}}_{\boldsymbol{Al}_{\boldsymbol{13}}\boldsymbol{Fe}_{\boldsymbol{4}}}$ [K]', fontdict=axfont)
    ax.set_ylabel(r'$\boldsymbol{\Delta}$ T$^{\boldsymbol{*}}_{\boldsymbol{Al}_{\boldsymbol{6}}\boldsymbol{Fe}}$ [K]', fontdict=axfont)
    ax.grid(True, linestyle='--', alpha=0.6)

    # # 4. Handle the legend
    # # Get the legend object created by seaborn
    # legend = ax.get_legend()
    
    # # Create descriptive labels for the legend
    # handles = legend.legend_handles
    # labels = [h.get_label() for h in handles]
    
    # # Replace short codes with full descriptions
    # descriptive_labels = [f"{code} : {description_map.get(code, 'N/A')}" for code in labels]
    
    # # Move the legend outside the plot area and update labels
    # ax.legend(handles, labels, 
    #           title="Configurations",
    #           bbox_to_anchor=(1.05, 1), 
    #           loc='upper left', 
    #           borderaxespad=0.)

    # 5. Final layout adjustment
    # This ensures that the external legend fits within the saved figure
    fig.tight_layout()
    
    # Save and show the plot
    # fig.savefig("Example11_Map.svg", bbox_inches='tight')
    plt.show()



if __name__ == "__main__":

    # Generate Undercooling grid

    # Create a meshgrid
    values = np.linspace(0, 50, 51)
    xx, yy = np.meshgrid(values, values)

    # Flatten the arrays and stack them as columns
    # We use 'F' to flatten column-by-column, which gives the desired order.
    grid = np.stack([xx.flatten('F'), yy.flatten('F')], axis=1)
    n,_=grid.shape
    print(f'Total {n} calculations')
    
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlFeSi_99Liu')

    # Parse database
    LiquidPhaseName= 'LIQUID'
    Database = eq.read_dat(datafile + '.dat')
    Condition = {
    'T': [1000]*n,  # Initial temperature, will be updated by the function
    'P': [1]*n,
    'Al': [0.975]*n,
    'Fe': [0.025]*n
    }
    CriticalUndercooling={
        'FCC_A1': [0.5]*n,
        'AL13FE4':list(grid[:,0]),
        'Al6Fe_AL6MN(s)': list(grid[:,1])
    }

    dT= 0.1
    UnitIn = ['C', 'atm', 'g']
    StartFromLiquidus = True
    res =eq.nucleoscheil_constituent_batch(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT=0.1,Unit=UnitIn)

    df = pl.DataFrame(res)
    ucs = pl.DataFrame({f'dT_{key}':CriticalUndercooling[key] for key in CriticalUndercooling.keys()})
    df= pl.concat([df,ucs],how="horizontal")
    df= df.to_pandas()
    print(df)


    #Plot
    label_map = {
        'FCC_A1 [g]': 'A',
        'AL13FE4 [g]': 'B',
        'Al6Fe_AL6MN(s) [g]': 'C',
        'Al6Fe_AL6MN(s)+FCC_A1 [g]': 'd',
        'AL13FE4+FCC_A1 [g]': 'e'
    }
    description_map = {v: k for k, v in label_map.items()}
    # Also add single configurations to the description map
    for k, v in label_map.items():
        description_map[v] = k

    processed_df = process_data(df, label_map)
    
    # Generate the plot
    plot_phase_map(processed_df, description_map)