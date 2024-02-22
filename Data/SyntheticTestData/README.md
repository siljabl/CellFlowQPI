## Synthetic data for testing and validating DIC and Continuity Metropolis

Generated with `generate_data.py` and `generate_field.py` in /CellFlowQPI/Code/SyntheticData.
Each folder contains:
- config.txt: configuration file
- _velocity_full.txt: full velocity field used to generate .tif files. This field covers a larger area than the final frames.
- /tif: folder containing test frames and corresponding velocity field
- /plots: folder containing plots of full data (including padding region), velocity field and intensity conservation
- /dic (remove this?): folder containing results from DIC
