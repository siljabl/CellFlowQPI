## Code to generate synthetic test data for the continuity metropolis algorithm

### Generate a velocity field
Run <code>python generate_field.py</code> to generate a velocity field. This also creates a new folder for the velocity field unless it already exists.

### Generate frames
Run <code>python generate_data.py {input image} {directory containing velocity field} [-n_frames] [-pad_width]</code> to generate cell layers moving according to a given velocity field. Input image should be a .tif and will be used as initial condition of the cell layer.


### TO DO:
- make possible to run outside this directory
- create virtual environment
- save config as dict
- is u, v saved correctly compared to .tif images? (I.e. different origin)
. improve resolution on RK?