#!/bin/bash

python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/uniform/ -u_max=8
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/uniform/ -u_max=16


python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/curl/ -u_max=1
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/curl/ -u_max=2
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/curl/ -u_max=4
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/curl/ -u_max=8
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/curl/ -u_max=16


python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_1/ -u_max=1
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_1/ -u_max=2
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_1/ -u_max=4
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_1/ -u_max=8
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_1/ -u_max=16

python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_2/ -u_max=1
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_2/ -u_max=2
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_2/ -u_max=4
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_2/ -u_max=8
python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif  ../../Data/SyntheticTestData/field_2/ -u_max=16
