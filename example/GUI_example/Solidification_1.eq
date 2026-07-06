{
  "format": "equilipy.eq",
  "modules": [
    {
      "amount_unit": "wt%",
      "batch_condition": {},
      "batch_condition_path": "",
      "batch_cpu_count": 0,
      "calculation_type": "single",
      "composition": [
        {
          "amount": "?100",
          "species": "Al"
        },
        {
          "amount": "5.5",
          "species": "Mg"
        },
        {
          "amount": "2.2",
          "species": "Si"
        }
      ],
      "delta_t": "5",
      "id": "session-1:solidification:1",
      "kind": "solidification",
      "liquid_phase": "LIQUID",
      "name": "Solidification#1",
      "nucleation_undercooling": {},
      "phase_names": [],
      "pressure": "1",
      "pressure_unit": "atm",
      "result": {},
      "result_columns": [],
      "results_table": {
        "headers": [],
        "rows": [],
        "visible": false
      },
      "script": {
        "condition": {
          "Al": 92.3,
          "Mg": 5.5,
          "P": 1.0,
          "Si": 2.2,
          "T": 2000.0
        },
        "database_path": "databases/AlCuMgSi_ORNL_FS83.dat",
        "delta_t": 5.0,
        "liquid_phase": "LIQUID",
        "mode": "single",
        "n_cpu": 1,
        "nucleation_undercooling": {},
        "phases": [],
        "solidification_model": "scheil",
        "start_from_liquidus": true,
        "units": [
          "C",
          "atm",
          "wt%"
        ]
      },
      "selected_phase_names": [],
      "solidification_model": "scheil",
      "start_from_liquidus": true,
      "temperature": "2000",
      "temperature_unit": "C",
      "transition_search": true
    }
  ],
  "scope": "module",
  "session": {
    "amount_unit": "wt%",
    "databases": [
      {
        "name": "AlCuMgSi_ORNL_FS83.dat",
        "original_path": "/Users/69e/Documents/Projects/2026VTO_Equilipy/example/GUI_example/databases/AlCuMgSi_ORNL_FS83.dat",
        "path": "databases/AlCuMgSi_ORNL_FS83.dat",
        "selected": true
      }
    ],
    "default_phase_names": [],
    "id": "session-1",
    "name": "Session#1",
    "pressure_unit": "atm",
    "result_column_defaults": {},
    "temperature_unit": "C"
  },
  "version": 1
}