#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import polars as pl

import equilipy as eq

if __name__ == "__main__":
    # Parse database
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, "database", "AlFeSi_99Liu_v84")
    DB = eq.read_dat(datafile + ".dat")

    # Set input data
    NPT = {
        "T": 1000,  # Initial temperature, will be updated by the function
        "P": 1,
        "Al": 0.975,
        "Fe": 0.025,
    }

    TargetPhase = "LIQUID"

    # Calculate Scheil cooling
    res = eq.scheil_cooling(TargetPhase, DB, NPT, delta_T=1, unit=["K", "atm", "g"])

    segments = res.scheil_constituent_segments(x="fs_w", y="T", y_unit="input")

    fig, ax = plt.subplots(figsize=(5, 4.5))
    titlefont = {"size": 18, "family": "Arial", "color": "#000000", "weight": "bold"}
    axfont = {"size": 16, "family": "Arial", "color": "#000000", "weight": "bold"}
    lfont = {
        "size": 14,
        "family": "Arial",
    }
    ax.set_title("Scheil prediction Al-2.5wt%Fe", fontdict=titlefont)
    ax.set_xlim(0, 0.16)
    ax.set_ylim(900, 1000)
    ax.set_ylabel(r"Temperature [K]", fontdict=axfont)
    ax.set_xlabel(r"Fraction of all solids [g/g]", fontdict=axfont)

    colors = ["k", "#4169e1", "#ff7f50"]
    for i, segment in enumerate(segments):
        ax.plot(
            segment.x,
            segment.y,
            "-",
            lw=2,
            color=colors[i % len(colors)],
            label=segment.label,
        )
    ax.plot([0.018298033, 0.03], [927.194274902343, 950], "k-", lw=0.5)
    ax.text(0.03, 950, "1.83 wt%", fontdict=lfont)

    ax.legend(prop=lfont)
    fig.tight_layout()
    print("Scheil constituents")
    print(pl.DataFrame(res.scheil_constituents))
    plt.show()
    #
    # print(df)
