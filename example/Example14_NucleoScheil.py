import os
import equilipy as eq

xfont = {"size": 10, "family": "Arial", "color": "#000000", "weight": "bold"}
yfont = {"size": 16, "family": "Arial", "color": "#000000", "weight": "bold"}
if __name__ == "__main__":
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, "database", "AlCuMgSi_ORNL_FS83")

    # Parse database
    liquid_phase_name = "LIQUID"
    database = eq.read_dat(datafile + ".dat", factsage_8_plus=True)
    condition = dict({"T": 900, "P": 1, "Al": 0.923, "Mg": 0.055, "Si": 0.022})
    critical_undercooling = {
        "FCC_A1": 0.5,
        "DIAMOND_A4": 1,
        "Mg2Si_MG2SI(s)": 2,
        "HCP_A3": 0,
    }

    res_nucleoscheil = eq.nucleoscheil_cooling(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        delta_T=0.1,
        unit=["C", "atm", "g"],
    )
    print(res_nucleoscheil.scheil_constituents)
