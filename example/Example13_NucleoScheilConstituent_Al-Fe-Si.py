import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import polars as pl
import equilipy as eq
xfont= {
    'size': 10,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
}
yfont= {
    'size': 16,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
}
if __name__ == "__main__":

    datafile = '../database/AlFeSi_99Liu'

    # Parse database
    LiquidPhaseName= 'LIQUID'
    Database = eq.read_dat(datafile + '.dat')
    Condition = {
        'T': 1000,  # Initial temperature, will be updated by the function
        'P': 1,
        'Al': 0.91,
        'Fe': 0.025,
        'Si': 0.065
    }
    CriticalUndercooling={
        'FCC_A1': 0.5,
        'ALFESI_ALPHA': 8,
        'ALFESI_BETA': 2, # caused by slight nucleation undercooling of FCC_A1
        'DIAMOND_A4': 1

    }
    res_scheil= res = eq.scheil_cooling(LiquidPhaseName,Database,Condition,dT=0.1,Unit=['C', 'atm', 'g']).ScheilConstituents

    res_nucleoscheil =eq.nucleoscheil_cooling(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT=0.1,Unit=['C', 'atm', 'g']).ScheilConstituents
    df = pl.concat([pl.DataFrame(res_scheil),pl.DataFrame(res_nucleoscheil)], how="diagonal")

    header = [constituent for constituent in df.columns if ' [g]' in constituent]
    header = header[3:]
    packages = ['Scheil', 'NucleoScheil']
    tasks = [x.split(' [g]')[0] for x in header]
        
    performance_data = np.zeros((len(header),2))

    for i, x in enumerate(header):
        try: 
            performance_data[i,0]=res_scheil[x]*100
        except KeyError:
            performance_data[i,0]=0
        try: 
            performance_data[i,1]=res_nucleoscheil[x]*100
        except KeyError:
            performance_data[i,1]=0
    # Replace these values with your actual performance data


    # Automatically calculate bar width based on the number of packages
    num_packages = len(packages)
    bar_width = 0.8 / num_packages  # 0.8 ensures some gap between grouped bars
    bar_positions = np.arange(len(tasks))

    # Colors and hatches
    colori = ['#ff7f50', '#4169e1']
    # colori = [ '#8060ff', 'r','b','g']
    hatches = ['xxx', 'oo', '///', '**', '\\\\\\', '...']



    # Create figure with GridSpec for precise control
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.03)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Colors and hatches
    colori = ['#ff7f50', '#4169e1']
    # colori = [ '#8060ff', 'r','b','g']
    hatches = ['xxx', 'oo', '///', '**', '\\\\\\', '...']

    # Plot the same data on both axes
    for ax in [ax1, ax2]:
        # Adding grid and labels
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')
        
        # Plotting bars dynamically
        for i in range(num_packages):
            name = packages[i]
            bars = ax.bar(bar_positions + i * bar_width, performance_data[:, i],
                        width=bar_width, label=packages[i] if ax == ax1 else "", 
                        edgecolor='k', color=colori[i])
            
            # Add value labels on top of bars - show on appropriate subplot based on value
            for j, bar in enumerate(bars):
                height = bar.get_height()
                # Show labels on the subplot that contains the bar
                if (ax == ax1 and height >= 30) or (ax == ax2 and height < 30):
                    # Position label above the bar with some padding
                    label_y = height + (0.1 if ax == ax1 else 0.1)
                    ax.text(bar.get_x() + bar.get_width()/2.+0.03, label_y,
                        f'{height:.0f}%', ha='center', va='bottom', 
                        fontsize=9, fontweight='bold')

    # Set different y-limits for broken axis
    ax1.set_ylim(30, 60)  # Top part (30-60)
    ax2.set_ylim(0, 9)   # Bottom part (0-10)
    ax2.set_yticks(np.arange(0, 9, 2)) # Ticks for bottom axis

    # Hide the spines between ax1 and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # Hide top x-axis ticks
    ax2.xaxis.tick_bottom()

    # Add break marks
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    # Set labels and ticks
    # Add a big axis, hide frame, and set the y-label
    ax_label = fig.add_subplot(111, zorder=-1) # zorder=-1 to put it behind other plots
    ax_label.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    ax_label.spines['top'].set_color('none')
    ax_label.spines['bottom'].set_color('none')
    ax_label.spines['left'].set_color('none')
    ax_label.spines['right'].set_color('none')
    ax_label.set_ylabel('Microstructural constituents, wt%', fontdict=yfont)


    # Center xticks in the middle of the group
    ax2.set_xticks(bar_positions + bar_width * (num_packages - 1) / 2)
    ax2.set_xticklabels(tasks, fontdict=xfont, rotation=10)

    # Adding legend (only to top subplot)
    ax1.legend(ncol=1, loc="upper left", prop={'size': 12})

    # Save and show
    fig.tight_layout()

    # Save and display
    # Save and display
    fig.savefig('ScheilConstituents.svg', transparent=False, bbox_inches='tight')
    plt.show()
