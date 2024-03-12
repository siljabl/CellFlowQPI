#include <string.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <math.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/non_central_chi_squared_distribution.hpp>
#include <omp.h>

void metropolis_hastings(double epsinit, double epsfinal, double epsdiff, int nosteps, double alpha, double temp, double phi1, double phi2, double gamma, double gamma2, double rhoinit[], double rhodiff[], double u[], double v[], double uc[], double vc[]);
