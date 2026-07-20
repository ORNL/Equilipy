import os

import matplotlib.pyplot as plt
import polars as pl
import equilipy as eq

if __name__ == "__main__":
    pl.Config.set_tbl_cols(10)
    pl.Config.set_tbl_rows(100)

    # datafile = '../database/AlFe_PanAl24'
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, "database", "AlFeSi_99Liu_v84")

    # Parse database
    liquid_phase_name = "LIQUID"
    database = eq.read_dat(datafile + ".dat")
    condition = {
        "T": 1000,  # Initial temperature, will be updated by the function
        "P": 1,
        "Al": 0.972,
        "Fe": 0.028,
    }
    critical_undercooling = {"FCC_A1": 0.5, "AL13FE4": 20, "Al6Fe_AL6MN(s)": 6}
    UnitIn = ["C", "atm", "g"]
    res = eq.nucleoscheil_cooling(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        delta_T=0.1,
        unit=UnitIn,
    )

    segments = res.scheil_constituent_segments(x="fs_w", y="T", y_unit="input")

    fig, ax = plt.subplots(figsize=(5, 4.5))
    titlefont = {"size": 18, "family": "Arial", "color": "#000000", "weight": "bold"}
    axfont = {"size": 16, "family": "Arial", "color": "#000000", "weight": "bold"}
    lfont = {
        "size": 14,
        "family": "Arial",
    }
    ax.set_title("Nucleoscheil prediction Al-2.8wt%Fe", fontdict=titlefont)
    ax.set_xlim(0, 0.16)
    # ax.set_ylim(920, 980)
    ax.set_ylabel(r"Temperature [C]", fontdict=axfont)
    ax.set_xlabel(r"Fraction of all solids [g/g]", fontdict=axfont)

    colors = ["#999999", "#4169e1", "#2e8b57", "#ff7f50"]
    for i, segment in enumerate(segments):
        ax.plot(
            segment.x,
            segment.y,
            "-",
            lw=2,
            color=colors[i % len(colors)],
            label=segment.label,
        )

    ax.legend(prop=lfont)
    fig.tight_layout()
    print("Nucleoscheil constituents")
    print(pl.DataFrame(res.scheil_constituents))
    plt.show()
