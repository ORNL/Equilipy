import numpy as np, pandas as pd, matplotlib.pyplot as plt

import polars as pl
import equilipy as eq
import os

if __name__ == "__main__":

    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlFeSi_99Liu')

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
    res_scheil= res = eq.scheil_cooling(LiquidPhaseName,Database,Condition,dT=0.1,Unit=['C', 'atm', 'g'])
    dict_scheil= res_scheil.to_dict()
    data_scheil_x = np.array(dict_scheil['LIQUID_Endmembers_Xmass_Si'])*100
    data_scheil_y = np.array(dict_scheil['LIQUID_Endmembers_Xmass_Fe'])*100
    data_scheil = np.asarray(list(zip(data_scheil_x,data_scheil_y)))

    res_nucleoscheil =eq.nucleoscheil_cooling(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT=0.1,Unit=['C', 'atm', 'g'])
    dict_nucleoscheil=res_nucleoscheil.to_dict()
    data_nucleoscheil_x = np.array(dict_nucleoscheil['LIQUID_Endmembers_Xmass_Si'])*100
    data_nucleoscheil_y = np.array(dict_nucleoscheil['LIQUID_Endmembers_Xmass_Fe'])*100
    data_nucleoscheil = np.asarray(list(zip(data_nucleoscheil_x,data_nucleoscheil_y)))

    axfont= {
        'size': 16,
        'family': 'Arial',
        'color' : '#000000',
        'weight': 'bold'
    }
    txtfont= {
        'size': 14,
        'family': 'Arial',
        'color' : '#000000',
        'ha': 'center',
        'va': 'center'
    }
    titlefont= {
    'size': 18,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
    }
    

    df=pd.read_csv(os.path.join(fpath, 'Projection_AlFeSi.csv'))
    header= df.keys()
    PDs= [h for h in header if '99Liu_' in h]
    # NucleoScheil = [h for h in header if 'NucleoScheil_' in h]
    # Scheil = [h for h in header if 'Scheil1_' in h]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4.5))
    
    # Plot solidification path
    ax.set_title('Solidification path',fontdict=titlefont)
    ax.set_xlabel('Si content [wt.%]',fontdict=axfont)
    ax.set_ylabel('Fe content [wt.%]',fontdict=axfont)
    ax.set_xlim(4, 14)
    ax.set_ylim(0.5, 3)

    # Plot Pandat data
    for i in range(int(len(PDs)/2)):
        data = df[[PDs[2*i],PDs[2*i+1]]].dropna().to_numpy()
        ax.plot(data[:,0],data[:,1], color='k', linestyle='-', linewidth=2)


    plt.legend(loc='upper left', fontsize=12, frameon=False)

    ax.plot(data_nucleoscheil[:,0],data_nucleoscheil[:,1], color='#4169e1', linestyle='-', linewidth=3)
    ax.plot(data_scheil[:,0],data_scheil[:,1], color='#ff7f50', linestyle='-', linewidth=3)

    ax.plot([6.5],[2.5],'o',color='k',ms=10,lw=3)
    ax.text(4.3,2.7,r'Al$_{13}$Fe$_4$',fontdict=txtfont, rotation=270)
    ax.text(7,2.9,r'$\alpha$-AlFeSi',fontdict=txtfont)
    ax.text(12,2.9,r'$\beta$-AlFeSi',fontdict=txtfont)
    ax.text(13.5,0.6,r'Si',fontdict=txtfont)
    ax.text(5,0.6,r'FCC_Al',fontdict=txtfont)
    fig.tight_layout()

    plt.show()