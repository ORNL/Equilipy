import os

import polars as pl

import _local_import  # noqa: F401

import equilipy as eq


if __name__ == "__main__":
    pl.Config.set_tbl_cols(40)
    pl.Config.set_tbl_rows(40)
    pl.Config.set_fmt_str_lengths(80)

    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, "database", "AlFeSi_99Liu_v84")

    liquid_phase_name = "LIQUID"
    database = eq.read_dat(datafile + ".dat")
    condition = {
        "T": 1000,
        "P": 1,
        "Al": 0.91,
        "Fe": 0.025,
        "Si": 0.065,
    }
    critical_undercooling = {
        "FCC_A1": 0.5,
        "ALFESI_ALPHA": 8,
        "ALFESI_BETA": 2,
        "DIAMOND_A4": 1,
    }

    result = eq.nucleoscheil_cooling(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        delta_T=0.1,
        unit=["C", "atm", "g"],
    )

    result_table = result.to_table()
    result_df = result_table.to_polars()
    constituents_df = pl.DataFrame(result.scheil_constituents)

    result_csv = os.path.join(fpath, "Result_Ex15_NucleoScheil_AlFeSi.csv")
    constituents_csv = os.path.join(
        fpath,
        "ScheilConstituents_Ex15_NucleoScheil_AlFeSi.csv",
    )
    result_table.to_csv(result_csv)
    constituents_df.write_csv(constituents_csv)

    preview_columns = [
        column
        for column in [
            "T [C]",
            "T [K]",
            "label",
            "fl_w",
            "fs_w",
            "Q [J]",
            "LIQUID_amount_w [g]",
            "FCC_A1_amount_w [g]",
            "ALFESI_BETA_amount_w [g]",
            "DIAMOND_A4_amount_w [g]",
        ]
        if column in result_df.columns
    ]
    preview_df = result_df.select(preview_columns)

    print("NucleoScheil result table")
    print(f"Rows: {result_df.height}, Columns: {result_df.width}")
    print("\nAvailable columns:")
    for column in result_table.available_columns():
        print(f"  - {column}")
    print("\nFocused first rows:")
    print(preview_df.head(10))
    print("\nFocused last rows:")
    print(preview_df.tail(10))

    print("\nNucleoScheil constituents")
    print(constituents_df)

    print("\nSaved files")
    print(result_csv)
    print(constituents_csv)
