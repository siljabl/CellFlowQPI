#!/bin/bash

python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/none/
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/uniform/
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/curl/
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/field_1/
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/field_2/
