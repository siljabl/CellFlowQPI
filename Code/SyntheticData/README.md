## Code to generate synthetic test data for the continuity metropolis algorithm

### Generate a velocity field
Run <code>python generate_field.py</code> to generate a velocity field. This also creates a new folder for the velocity field unless it already exists. The velocity field is first created as a function in src/velocity_fields.py and must be chosen as function in the script. generate_field.py creates a grid and hands it to the function in order to plot and save a given field.

### Generate frames
Run <code>python generate_data.py {input image} {directory containing velocity field} [-n_frames] [-pad_width] [-u_max]</code> to generate cell layers moving according to a given velocity field. Input image should be a .tif and will be used as initial condition of the cell layer.


### TO DO:
- make possible to run outside this directory
- create virtual environment
- save config as dict
- improve resolution on RK?
- remove normalization of velocity field? Since normalizing in generate_data.py, full velocity field and tif velocity is not equal...
- nb. tif_size not correct
