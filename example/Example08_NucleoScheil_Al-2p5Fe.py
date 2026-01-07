
import polars as pl
import equilipy as eq
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    pl.Config.set_tbl_cols(10)
    pl.Config.set_tbl_rows(100)

    # datafile = '../database/AlFe_PanAl24'
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlFeSi_99Liu')

    # Parse database
    LiquidPhaseName= 'LIQUID'
    Database = eq.read_dat(datafile + '.dat')
    Condition = {
    'T': 1000,  # Initial temperature, will be updated by the function
    'P': 1,
    'Al': 0.975,
    'Fe': 0.025
    }
    CriticalUndercooling={
        'FCC_A1': 0.5,
        'AL13FE4':13.5
    }
    UnitIn = ['C', 'atm', 'g']
    res = eq.nucleoscheil_cooling(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT=0.1,Unit=UnitIn)
    


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
    ax.set_title('NucleoScheil prediction Al-2.5wt%Fe',fontdict=titlefont)
    ax.set_xlim(0,0.16)
    ax.set_ylim(920,980)
    ax.set_ylabel(r'Temperature [K]',fontdict=axfont )
    ax.set_xlabel(r'Fraction of all solids [g/g]',fontdict=axfont)

    colors = ['#999999','#4169e1','#2e8b57','#ff7f50']
    for i, label in enumerate(labels):
        fs_g,T_K = dataset[i]
        ax.plot(fs_g,T_K,'-',lw=2,color=colors[i],label=label)
    
    ax.plot([0.018511637,0.03],[926.746463979682,922],'k-',lw=0.5)
    ax.text(0.03,922,'1.85 wt%',fontdict=lfont)

    ax.plot([0.137730888046083,0.133],[926.246463979682,925],'k-',lw=0.5)
    ax.text(0.12,922,'13.77 wt%',fontdict=lfont)

    ax.legend(prop=lfont)
    fig.tight_layout()
    print('NucleoScheil constituents')
    print(pl.DataFrame(res.ScheilConstituents))
    plt.show()
    