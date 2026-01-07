#!/usr/bin/env python3
import polars as pl
import equilipy as eq
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    #Parse database
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlFeSi_99Liu')
    DB=eq.read_dat(datafile+'.dat')

    # Set input data
    NTP = {
    'T': 1000,  # Initial temperature, will be updated by the function
    'P': 1,
    'Al': 0.975,
    'Fe': 0.025
    }

    TargetPhase='LIQUID'

    # Calculate Scheil cooling
    res=eq.scheil_cooling(TargetPhase,DB,NTP,dT=1,Unit=['C','atm','g'])

    n = len(res.PhaseLabel)
    dataset = []
    
    fs_g = [res.fs_mass[0]]
    T_K = [res.T[0]]
    labels = [res.PhaseLabel[0]]
    for i in range(1,n):
        label = res.PhaseLabel[i]
        if label != labels[-1]:
            
            labels.append(label)
            dataset.append([fs_g,T_K])
            #Initialize dataset
            fs_g = [fs_g[-1]]
            T_K = [T_K[-1]]
            
        fs_g.append(res.fs_mass[i])
        T_K.append(res.T[i])
    dataset.append([fs_g,T_K])

    fig, ax = plt.subplots(figsize=(5,4.5))
    titlefont= {
    'size': 18,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
    }
    axfont= {
    'size': 16,
    'family': 'Arial',
    'color' : '#000000',
    'weight': 'bold'
    }
    lfont = {
        'size': 14,
        'family': 'Arial',
    }
    ax.set_title('Scheil prediction Al-2.5wt%Fe',fontdict=titlefont)
    ax.set_xlim(0,0.16)
    ax.set_ylim(900,1000)
    ax.set_ylabel(r'Temperature [K]',fontdict=axfont )
    ax.set_xlabel(r'Fraction of all solids [g/g]',fontdict=axfont)

    colors = ['k','#4169e1','#ff7f50']
    for i, label in enumerate(labels):
        
        fs_g,T_K = dataset[i]
        ax.plot(fs_g,T_K,'-',lw=2,color=colors[i],label=label)
    ax.plot([0.018298033,0.03],[927.194274902343,950],'k-',lw=0.5)
    ax.text(0.03,950,'1.83 wt%',fontdict=lfont)
    
    ax.legend(prop=lfont)
    fig.tight_layout()
    print('Scheil constituents')
    print(pl.DataFrame(res.ScheilConstituents))
    plt.show()
    # 
    # print(df)
