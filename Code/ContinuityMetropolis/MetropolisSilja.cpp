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

//#define NX 196 // 1/4 of image size, change this!
#define PI 3.14159265359

using namespace std;

int main(int argc, char** argv) {

	int res = 4;	// number of input pixels per output pixels
	int image_x_res = atoi(argv[5]); // read from image or from python
	int image_y_res = atoi(argv[6]); // read from image or from python

	int nx = (int) image_x_res / res;
	int ny = (int) image_y_res / res;
	cout << "resolution: (" << nx << ", " << ny << ")" << endl;

	double *imagerhoinit, *imagerhofinal, *rhoinit, *rhofinal, *rhodiff, *u, *v, *x1, *x2, *imageuc, *imagevc, *uc, *vc;

// ***********************************************
//	Memory allocation
// ***********************************************
	
	imagerhoinit  = new double[image_x_res * image_y_res];
	imagerhofinal = new double[image_x_res * image_y_res];
	
	rhoinit  = new double[nx * ny];
	rhofinal = new double[nx * ny];
	rhodiff  = new double[nx * ny];
		
	u = new double[nx * ny];
	v = new double[nx * ny];
	
	x1 = new double[nx * ny];
	x2 = new double[nx * ny];
	
	imageuc = new double[image_x_res * image_y_res];
	imagevc = new double[image_x_res * image_y_res];
	
	uc = new double[nx * ny];
	vc = new double[nx * ny];


// ***********************************************
//	Parameter setting
// ***********************************************

	double epsinit  = 0.0;
	double epsfinal = 0.0;
	double epsdiff  = 0.0;
		
	int seed1 = atoi(argv[7]);
	int nosteps = atoi(argv[8]);
	
	double alpha = atof(argv[9]);
	double temp  = atof(argv[10]);
	double phi1  = atof(argv[11]);
	double phi2  = atof(argv[12]);
	double gamma  = atof(argv[13]);
	double gamma2 = atof(argv[14]);
	
// ***********************************************
//	Noise initialisation
// ***********************************************

	// these are 1D, how to treat nx if nx != ny
	std::default_random_engine generator;
	
 	std::normal_distribution<double> distribution(0.0,0.1);
 	std::uniform_real_distribution<double> distribution2(-1.0,1.0);
 	std::uniform_real_distribution<double> distribution3(0.0,1.0);
 	std::uniform_int_distribution<> distribution4(0, nx-1);
 	
 	std::uniform_int_distribution<> distributionx(2, nx-2);
 	std::uniform_int_distribution<> distributiony(2, ny-2);
 	
 	generator.seed(seed1);
 	
 	// delete the below and use: 
 	// std::normal_distribution<double> dist(0.0,0.05);
 	// If you can't install boost
 	
	boost::mt19937 gen;
	boost::random::normal_distribution<> dist(0.0,0.05);
	gen.seed(seed1);
	dist.reset();

// ***********************************************
//	Initialisation of values
// ***********************************************

	double xmod, ymod, radterm, x1temp, x2temp;
	

	// --- fixed initialisation of the fields
	for (int i = 0; i < nx; i++) {		
		for (int j = 0; j < ny; j++) {
		
			rhoinit[i * ny + j]  = 0.0;
			rhofinal[i * ny + j] = 0.0;
			rhodiff[i * ny + j]  = 0.0;
			
			u[i * ny + j] = 0.0;
			v[i * ny + j] = 0.0;
			
			uc[i * ny + j] = 0.0;
			vc[i * ny + j] = 0.0;
			
			xmod = (i * 1.0);	// remove this?
			ymod = (j * 1.0);	// remove this?
			
			x1[i * ny + j] = xmod;
			x2[i * ny + j] = ymod;
		}
	}


// ***********************************************
//	Reading of Input
// ***********************************************	
		
	string inImage1 = argv[1];
	string inImage2 = argv[2];
	string inDICx = argv[3];
	string inDICy = argv[4];

    ifstream Image1, Image2, xDIC, yDIC;

    Image1.open(inImage1.c_str());
    Image2.open(inImage2.c_str());
    xDIC.open(inDICx.c_str());
    yDIC.open(inDICy.c_str());

	for (int n = 0; n < image_x_res; n++) {
		for (int m = 0; m < image_y_res; m++) {
			Image1 >>  imagerhoinit[n * image_x_res + m];
			Image2 >> imagerhofinal[n * image_x_res + m];

			xDIC >> imageuc[n * image_x_res + m];
			yDIC >> imagevc[n * image_x_res + m];
		}		
	}
 	// close input files
    Image1.close();
    Image2.close();
	xDIC.close();
	yDIC.close();

// ***********************************************
//	Interpolate to used grid
// ***********************************************

	int usedn, usedm;
	
	for (int i = 0; i < nx; i++) {		
		for (int j = 0; j < ny; j++) {
			for (int n = 0; n < res; n++) {
				for (int m = 0; m < res; m++) {
					usedn = i * res + n;
					usedm = j * res + m;
					rhoinit[i * ny + j]  =  rhoinit[i * ny + j] +  imagerhoinit[usedn * image_y_res + usedm];
					rhofinal[i * ny + j] = rhofinal[i * ny + j] + imagerhofinal[usedn * image_y_res + usedm];
					
					uc[i * ny + j] = uc[i * ny + j] + gamma2 * imageuc[usedn * image_y_res + usedm];
					vc[i * ny + j] = vc[i * ny + j] + gamma2 * imagevc[usedn * image_y_res + usedm];
				}
			}
			
		}
	}
	
	for (int i = 0; i < nx; i++) {		
		for (int j = 0; j < ny; j++) {
			int idx = i * ny + j;
			
			rhodiff[idx] = rhofinal[idx] - rhoinit[idx];
			rhoinit[idx] = 0.5 * rhofinal[idx] + 0.5 * rhoinit[idx];
			
			u[idx] = uc[idx];
			v[idx] = vc[idx];
		}
	}
	

// ***********************************************
//	Metropolis-Hastings
// ***********************************************

	int chosenindex = 0;
	double probofmove = 0.0;
	
	int i = 0;
	int j = 0;
	
	int idx0, idx1, idx2, idx3, idx4;
	int index1i, index1j, index2i, index2j, index3i, index3j, index4i, index4j;
	double initJd, initJr, initJc, potJd, potJr, potJc;
	
	double potJc1, potJc2, initJc1, initJc2;
	
	double initJd1, initJd2, initJd3, initJd4, initJd0;
	double potJd1, potJd2, potJd3, potJd4, potJd0;
	double initJd_tmp, potJd_tmp;
	
	double init1dux, init1dvx, init2dux, init2dvx, init3dux, init3dvx, init4dux, init4dvx;
	double init1duy, init1dvy, init2duy, init2dvy, init3duy, init3dvy, init4duy, init4dvy;
	double pot1dux, pot1dvx, pot2dux, pot2dvx, pot3dux, pot3dvx, pot4dux, pot4dvx;
	double pot1duy, pot1dvy, pot2duy, pot2dvy, pot3duy, pot3dvy, pot4duy, pot4dvy;
	double uinc, vinc;
	
	double trueenergy;

	int numchanges = 0;
	
	for (int s = 0; s < nosteps; s++) {
	
		// only check energy of changed values
		
		i = distributionx(generator);
		j = distributiony(generator);

		idx0 = i * ny + j;
		
		index1i = i - 1;
		index1j = j;
		idx1 = index1i * ny + index1j;
		
		index2i = i + 1;
		index2j = j;
		idx2 = index2i * ny + index2j;
		
		index3i = i;
		index3j = j + 1;
		idx3 = index3i * ny + index3j;
		
		index4i = i;
		index4j = j - 1;
		idx4 = index4i * ny + index4j;
		
		
		initJd = 0;
		int indexi[5] = {i, index1i, index2i, index3i, index4i};
		int indexj[5] = {j, index1j, index2j, index3j, index4j};
		
		for (int jdi = 0; jdi < 5; jdi++) {
			int index = indexi[jdi] * ny + indexj[jdi];

			initJd_tmp = abs(rhodiff[index] + 
							 phi1 * rhoinit[index] * (u[index + ny] - u[index - ny]) + 
							 phi1 * rhoinit[index] * (v[index + 1]  - v[index - 1])  +
							 phi2 * u[index] * (rhoinit[index + ny] - rhoinit[index - ny]) + 
							 phi2 * v[index] * (rhoinit[index + 1]  - rhoinit[index - 1]));

			initJd += initJd_tmp;
		}
		uinc = dist(gen);
		vinc = dist(gen);

		potJd = initJd + phi1 * uinc * (rhoinit[idx1] - rhoinit[idx2])
					   - phi1 * vinc * (rhoinit[idx3] - rhoinit[idx4])
					   + phi2 * uinc * (rhoinit[idx0 + ny] - rhoinit[idx0 - ny])
					   + phi2 * vinc * (rhoinit[idx0 + 1]  - rhoinit[idx0 - 1]);


			
		// IDK why, but doing this with the 9-stencil or whatever, gives terrible results, OR: the results are preferentially diagonal... idk why
		
		init1dux = u[idx1 + ny] - u[idx1 - ny];
		init1dvx = v[idx1 + ny] - v[idx1 - ny];
		init1dvy = v[idx1 + 1]  - v[idx1 - 1];
		init1duy = u[idx1 + 1]  - u[idx1 - 1];
		
		init2dux = u[idx2 + ny] - u[idx2 - ny];
		init2dvx = v[idx2 + ny] - v[idx2 - ny];
		init2dvy = v[idx2 + 1]  - v[idx2 - 1];
		init2duy = u[idx2 + 1]  - u[idx2 - 1];
		
		init3dux = u[idx3 + ny] - u[idx3 - ny];
		init3dvx = v[idx3 + ny] - v[idx3 - ny];
		init3dvy = v[idx3 + 1]  - v[idx3 - 1];
		init3duy = u[idx3 + 1]  - u[idx3 - 1];
		
		init4duy = u[idx4 + ny] - u[idx4 - ny];
		init4dvx = v[idx4 + ny] - v[idx4 - ny];
		init4dvy = v[idx4 + 1]  - v[idx4 - 1];
		init4duy = u[idx4 + 1]  - u[idx4 - 1];
		
		initJr = init1dux * init1dux + init1duy * init1duy + init1dvx * init1dvx + init1dvy * init1dvy +
				 init2dux * init2dux + init2duy * init2duy + init2dvx * init2dvx + init2dvy * init2dvy +
				 init3dux * init3dux + init3duy * init3duy + init3dvx * init3dvx + init3dvy * init3dvy +
				 init4dux * init4dux + init4duy * init4duy + init4dvx * init4dvx + init4dvy * init4dvy;
		
	//	initJc = (u[i * ny + j] - uc[i * ny + j]) * (u[i * ny + j] - uc[i * ny + j]) 
	//		+ (v[i * ny + j] - vc[i * ny + j]) * (v[i * ny + j] - vc[i * ny + j]);
	
	//	initJc1 = gamma * exp((u[i * ny + j] - uc[i * ny + j]) * (u[i * ny + j] - uc[i * ny + j]));
	//	initJc2 = gamma * exp((v[i * ny + j] - vc[i * ny + j]) * (v[i * ny + j] - vc[i * ny + j]));
		initJc1 = gamma * (u[idx0] - uc[idx0]) * (u[idx0] - uc[idx0]);
		initJc2 = gamma * (v[idx0] - vc[idx0]) * (v[idx0] - vc[idx0]);
		initJc = initJc1 + initJc2;
		
		epsinit = initJd + alpha * initJr + initJc;
	//	epsinit = initJd + alpha * initJr + gamma * exp(initJc);
	//	epsinit = initJd + alpha * initJr;

			
		pot1dux = u[idx1 + ny] - u[idx1 - ny] + uinc;
		pot1dvx = v[idx1 + ny] - v[idx1 - ny] + vinc;
		pot1dvy = v[idx1 + 1] - v[idx1 - 1];
		pot1duy = u[idx1 + 1] - u[idx1 - 1];
		
		pot2dux = u[idx2 + ny] - u[idx2 - ny] - uinc;
		pot2dvx = v[idx2 + ny] - v[idx2 - ny] - vinc;
		pot2dvy = v[idx2 + 1] - v[idx2 - 1];
		pot2duy = u[idx2 + 1] - u[idx2 - 1];
		
		pot3dux = u[idx3 + ny] - u[idx3 - ny];
		pot3dvx = v[idx3 + ny] - v[idx3 - ny];
		pot3dvy = v[idx3 + 1] - v[idx3 - 1] - vinc;
		pot3duy = u[idx3 + 1] - u[idx3 - 1] - uinc;
		
		pot4dux = u[idx4 + ny] - u[idx4 - ny];
		pot4dvx = v[idx4 + ny] - v[idx4 - ny];
		pot4dvy = v[idx4 + 1] - v[idx4 - 1] + vinc;
		pot4duy = u[idx4 + 1] - u[idx4 - 1] + uinc;
		
		potJr = pot1dux * pot1dux + pot1duy * pot1duy + pot1dvx * pot1dvx + pot1dvy * pot1dvy +
				pot2dux * pot2dux + pot2duy * pot2duy + pot2dvx * pot2dvx + pot2dvy * pot2dvy +
				pot3dux * pot3dux + pot3duy * pot3duy + pot3dvx * pot3dvx + pot3dvy * pot3dvy +
				pot4dux * pot4dux + pot4duy * pot4duy + pot4dvx * pot4dvx + pot4dvy * pot4dvy;
		
	//	potJc = (u[i * ny + j] - uc[i * ny + j] + uinc) * (u[i * ny + j] - uc[i * ny + j] + uinc) 
	//		+ (v[i * ny + j] - vc[i * ny + j] + vinc) * (v[i * ny + j] - vc[i * ny + j] + vinc);
		
	//	potJc1 = gamma * exp((u[i * ny + j] - uc[i * ny + j] + uinc) * (u[i * ny + j] - uc[i * ny + j] + uinc));
	//	potJc2 = gamma * exp((v[i * ny + j] - vc[i * ny + j] + vinc) * (v[i * ny + j] - vc[i * ny + j] + vinc));
		potJc1 = gamma * (u[idx0] - uc[idx0] + uinc) * (u[idx0] - uc[idx0] + uinc);
		potJc2 = gamma * (v[idx0] - vc[idx0] + vinc) * (v[idx0] - vc[idx0] + vinc);
		potJc = potJc1 + potJc2;
		
		epsfinal = potJd + alpha * potJr + potJc;
	//	epsfinal = potJd + alpha * potJr + gamma * exp(potJc);
	//	epsfinal = potJd + alpha * potJr;
		
		// --- decide on swap --- 
		epsdiff = epsfinal - epsinit;
		epsdiff = temp * epsdiff;
		
//		cout << epsinit << "\t" << epsfinal << "\t" << epsdiff << endl;
		
		probofmove = exp(-epsdiff);
		
		if (epsdiff < 0.0) {
			u[idx0] = u[idx0] + uinc;
			v[idx0] = v[idx0] + vinc;
			numchanges = numchanges + 1;
		}
		else{
			if (probofmove > distribution3(generator)) {
				u[idx0] = u[idx0] + uinc;
				v[idx0] = v[idx0] + vinc;
				numchanges = numchanges + 1;
			}
		} 
		
	
	//

	/*
	if ((s % 10000) == 0) {				//Checking sum of energy
		trueenergy = 0.0;
	//	for (int m = 2; m < nx - 2; m++) {		
	//		for (int n = 2; n < nx - 2; n++) {
		for (int m = 3; m < nx - 3; m++) {		
			for (int n = 3; n < nx - 3; n++) {
			
				initJd = abs(rhodiff[m * nx + n] + 
					phi1 * rhoinit[m * nx + n] * (u[(m + 1) * nx + n] - u[(m - 1) * nx + n]) +
					phi1 * rhoinit[m * nx + n] * (v[m * nx + n + 1] - v[m * nx + n - 1]) +
					phi2 * (rhoinit[(m + 1) * nx + n] - rhoinit[(m - 1) * nx + n]) * u[m * nx + n] +
					phi2 * (rhoinit[m * nx + n + 1] - rhoinit[m * nx + n - 1]) * v[m * nx + n]);
				init1dux = u[(m + 1) * nx + n] - u[(m - 1) * nx + n];
				init1dvx = v[(m + 1) * nx + n] - v[(m - 1) * nx + n];
				init1dvy = v[m * nx + n + 1] - v[m * nx + n - 1];
				init1duy = u[m * nx + n + 1] - u[m * nx + n - 1];
				initJr = init1dux * init1dux + init1duy * init1duy + init1dvx * init1dvx + init1dvy * init1dvy;
				
			//	initJc = (u[m * nx + m] - uc[m * nx + n]) * (u[m * nx + n] - uc[m * nx + n]) 
			//	+ (v[m * nx + n] - vc[m * nx + n]) * (v[m * nx + n] - vc[m * nx + n]);
		
			//	initJc1 = gamma * exp((u[m * nx + n] - uc[m * nx + n]) * (u[m * nx + n] - uc[m * nx + n));
			//	initJc2 = gamma * exp((v[m * nx + n] - vc[m * nx + n]) * (v[m * nx + n] - vc[m * nx + n]));
				initJc1 = gamma * (u[m * nx + n] - uc[m * nx + n]) * (u[m * nx + n] - uc[m * nx + n]);
				initJc2 = gamma * (v[m * nx + n] - vc[m * nx + n]) * (v[m * nx + n] - vc[m * nx + n]);
				initJc = initJc1 + initJc2;
		
			//	trueenergy = trueenergy + initJd + alpha * initJr;
			
				trueenergy = trueenergy + initJd + alpha * initJr + initJc;
				
			//	trueenergy = trueenergy + initJd + alpha * initJr + gamma * initJc;
			}
		}
	
		cout << s << "\t" << trueenergy << endl;
	}  
		*/
	
	} // --- end of stepping ---
	
	
	/*
	for (int i = 0; i < 100; i++) {		
		for (int j = 0; j < 100; j++) {
			xmod = (i * 1.0) - 50;
			ymod = (j * 1.0) - 50;
			cout << xmod << "\t" << ymod << "\t" << rhoinit[i * 100 + j] << "\t" << rhofinal[i * 100 + j] << "\t" << rhodiff[i * 100 + j] << "\t" << j1[i * 100 + j] << "\t" << j2[i * 100 + j] << endl;
		}
	}
	*/
	
//	cout << "this many changes: " << numchanges << endl;

// ***********************************************
//	Output
// ***********************************************
	
	double Px, Py, pointx, pointy, sumweight, weight, pcdist, IntDensity;
	double bindist = 2.0; // same as out_res?

	int out_res = 4;
	int nx_out = (int) nx / out_res;
	int ny_out = (int) ny / out_res;

	ofstream xOutput("output/x_displacement.txt");
	ofstream yOutput("output/y_displacement.txt");
	ofstream iOutput("output/intensity.txt");

	// first we make an interpolation of the vector field
        for (int i = 0; i < nx_out; i++) {
        	for (int j = 0; j < ny_out; j++) {
        	
        		pointx = i * out_res;
        		pointy = j * out_res;
        		
        		Px = 0.0;
        		Py = 0.0;
        		
        		sumweight = 0.0;
        		weight = 1.0;
        		
        		IntDensity = 0.0;
        	
        		for (int n = 0; n < (nx * ny); n++){ 
        			pcdist = (x1[n] - pointx) * (x1[n] - pointx) + (x2[n] - pointy) * (x2[n] - pointy) + 1e-6;

        			if (pcdist < (bindist*bindist)){		// this 6 distance is hard coded, it shouldn't be
        		//		weight = sqrt(pcdist);       //doing these square roots is painful
        		//		weight = (bindist - weight)*(bindist - weight) / (bindist * bindist * pcdist);	//also doing these divisions
        				weight = 1.0;
        				Px = Px + weight * u[n];
        				Py = Py + weight * v[n];
        				sumweight = sumweight + weight;
        				
        				IntDensity = IntDensity + weight * rhoinit[n];
        				
        			}
        		}
        		
        		if(sumweight > 0.0) {		
        			Px = Px / sumweight;	
        			Py = Py / sumweight;
        			IntDensity = IntDensity / sumweight;
        		}
				xOutput << Px << " ";
				yOutput << Py << " ";
				iOutput << IntDensity << " ";
        	}
			xOutput << endl;
			yOutput << endl;
			iOutput << endl;
        }

	xOutput.close();
	yOutput.close();
	iOutput.close();
	


}	// --- end of program ---




