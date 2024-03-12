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

#define NX 256
#define PI 3.14159265359

// to compile: g++ -std=c++11 -O3 -fopenmp -o output.x MetropolisVer8.cpp

using namespace std;

int main(int argc, char** argv) {

	int res = 4;
	int image_x_res = 1024;
	int image_y_res = 1024;			

	int nx = NX;
	
	double *imagerhoinit, *imagerhofinal, *rhoinit, *rhofinal, *rhodiff, *partialIx, *partialIy, *partialIt, *u, *v, *x1, *x2, *imageuc, *imagevc, *uc, *vc;

// ***********************************************
//	Memory allocation
// ***********************************************
	
	imagerhoinit = new double[image_x_res * image_y_res];
	imagerhofinal = new double[image_x_res * image_y_res];
	
	rhoinit = new double[nx * nx];
	rhofinal = new double[nx * nx];
	rhodiff = new double[nx * nx];
		
	u = new double[nx * nx];
	v = new double[nx * nx];
	
	x1 = new double[nx * nx];
	x2 = new double[nx * nx];
	
	imageuc = new double[image_x_res * image_y_res];
	imagevc = new double[image_x_res * image_y_res];
	
	uc = new double[nx * nx];
	vc = new double[nx * nx];
	
// ***********************************************
//	Parameter setting
// ***********************************************

	double epsinit  = 0.0;
	double epsfinal = 0.0;
	double epsdiff  = 0.0;
		
	int seed1 = atoi(argv[1]);
	int nosteps = atoi(argv[2]);
	
	double alpha = atof(argv[3]);
	double temp = atof(argv[4]);
	double phi1 = atof(argv[5]);
	double phi2 = atof(argv[6]);
	double gamma = atof(argv[7]);
	double gamma2 = atof(argv[8]);
	
// ***********************************************
//	Noise initialisation
// ***********************************************

	std::default_random_engine generator;
	
 	std::normal_distribution<double> distribution(0.0,0.1);
 	std::uniform_real_distribution<double> distribution2(-1.0,1.0);
 	std::uniform_real_distribution<double> distribution3(0.0,1.0);
 	std::uniform_int_distribution<> distribution4(0, NX-1);
 	
 	std::uniform_int_distribution<> distributionx(2, NX-2);
 	std::uniform_int_distribution<> distributiony(2, NX-2);
 	
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
	for (int i = 0; i < NX; i++) {		
		for (int j = 0; j < NX; j++) {
		
			rhoinit[i * NX + j] = 0.0;
			rhofinal[i * NX + j] = 0.0;
			rhodiff[i * NX + j] = 0.0;
			
			u[i * NX + j] = 0.0;
			v[i * NX + j] = 0.0;
			
			uc[i * NX + j] = 0.0;
			vc[i * NX + j] = 0.0;
			
			xmod = (i * 1.0);
			ymod = (j * 1.0);
			
			x1[i * NX + j] = xmod;
			x2[i * NX + j] = ymod;
			
		}
	}

// ***********************************************
//	Reading of Input
// ***********************************************	

	int readXint, readYint;
		
	string inFileName = "image1.txt";
        ifstream inFile;
        inFile.open(inFileName.c_str());
	
	for (int n = 0; n < image_x_res; n++) {
		for (int m = 0; m < image_y_res; m++) {	
			inFile >> readXint >> readYint >> imagerhoinit[m * image_y_res + (image_y_res - 1 - n)];
		}		
	}
	
        inFile.close(); // close input file
        
        inFileName = "image2.txt";
        inFile.open(inFileName.c_str());
	
	for (int n = 0; n < image_x_res; n++) {
		for (int m = 0; m < image_y_res; m++) {	
			inFile >> readXint >> readYint >> imagerhofinal[m * image_y_res + (image_y_res - 1 - n)];
		}		
	}
	
        inFile.close(); // close input file
        
        inFileName = "dx1.txt";
        inFile.open(inFileName.c_str());
	
	for (int n = 0; n < image_x_res; n++) {
		for (int m = 0; m < image_y_res; m++) {	
			inFile >> readXint >> readYint >> imageuc[m * image_y_res + (image_y_res - 1 - n)];
		}		
	}
	
        inFile.close(); // close input file
        
        inFileName = "dy1.txt";
        inFile.open(inFileName.c_str());
	
	for (int n = 0; n < image_x_res; n++) {
		for (int m = 0; m < image_y_res; m++) {	
			inFile >> readXint >> readYint >> imagevc[m * image_y_res + (image_y_res - 1 - n)];
		}		
	}
	
        inFile.close(); // close input file


// ***********************************************
//	Interpolate to used grid
// ***********************************************

	int usedn, usedm;
	
	for (int i = 0; i < NX; i++) {		
		for (int j = 0; j < NX; j++) {
			for (int n = 0; n < res; n++) {
				for (int m = 0; m < res; m++) {			//honestly, this will probably only work if the image is square TODO fix this
					usedn = i * res + n;
					usedm = j * res + m;
					rhoinit[i * NX + j] = rhoinit[i * NX + j] + imagerhoinit[usedn * image_y_res + usedm];
					rhofinal[i * NX + j] = rhofinal[i * NX + j] + imagerhofinal[usedn * image_y_res + usedm];
					
					uc[i * NX + j] = uc[i * NX + j] + gamma2 * imageuc[usedn * image_y_res + usedm];
					vc[i * NX + j] = vc[i * NX + j] + gamma2 * imagevc[usedn * image_y_res + usedm];
				}
			}
			
		}
	}
	
	for (int i = 0; i < NX; i++) {		
		for (int j = 0; j < NX; j++) {
			rhodiff[i * NX + j] = rhofinal[i * NX + j] - rhoinit[i * NX + j];
			rhoinit[i * NX + j] = 0.5 * rhofinal[i * NX + j] + 0.5 * rhoinit[i * NX + j];;
			
			u[i * NX + j] = uc[i * NX + j];
			v[i * NX + j] = vc[i * NX + j];
		}
	}
	

// ***********************************************
//	Metropolis-Hastings
// ***********************************************

	int chosenindex = 0;
	double probofmove = 0.0;
	
	int i = 0;
	int j = 0;
	
	int index1i, index1j, index2i, index2j, index3i, index3j, index4i, index4j;
	double initJd, initJr, initJc, potJd, potJr, potJc;
	
	double potJc1, potJc2, initJc1, initJc2;
	
	double initJd1, initJd2, initJd3, initJd4, initJd0;
	double potJd1, potJd2, potJd3, potJd4, potJd0;
	
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
		
		index1i = i - 1;
		index1j = j;
		
		index2i = i + 1;
		index2j = j;
		
		index3i = i;
		index3j = j + 1;
		
		index4i = i;
		index4j = j - 1;
		
		initJd1 = abs(rhodiff[index1i * NX + index1j] + 
			phi1 * rhoinit[(index1i) * NX + index1j] * (u[(index1i + 1) * NX + index1j] - u[(index1i - 1) * NX + index1j]) +
			phi1 * rhoinit[index1i * NX + index1j] * (v[index1i * NX + index1j + 1] - v[index1i * NX + index1j - 1]) +
			phi2 * (rhoinit[(index1i + 1) * NX + index1j] - rhoinit[(index1i - 1) * NX + index1j]) * u[index1i * NX + index1j] +
			phi2 * (rhoinit[index1i * NX + index1j + 1] - rhoinit[index1i * NX + index1j - 1]) * v[index1i * NX + index1j]
			);
		initJd2 = abs(rhodiff[index2i * NX + index2j] + 
			phi1 * rhoinit[(index2i) * NX + index2j] * (u[(index2i + 1) * NX + index2j] - u[(index2i - 1) * NX + index2j]) +
			phi1 * rhoinit[index2i * NX + index2j] * (v[index2i * NX + index2j + 1] - v[index2i * NX + index2j - 1]) +
			phi2 * (rhoinit[(index2i + 1) * NX + index2j] - rhoinit[(index2i - 1) * NX + index2j]) * u[index2i * NX + index2j] +
			phi2 * (rhoinit[index2i * NX + index2j + 1] - rhoinit[index2i * NX + index2j - 1]) * v[index2i * NX + index2j]
			);	
		initJd3 = abs(rhodiff[index3i * NX + index3j] + 
			phi1 * rhoinit[(index3i) * NX + index3j] * (u[(index3i + 1) * NX + index3j] - u[(index3i - 1) * NX + index3j]) +
			phi1 * rhoinit[index3i * NX + index3j] * (v[index3i * NX + index3j + 1] - v[index3i * NX + index3j - 1]) + 
			phi2 * (rhoinit[(index3i + 1) * NX + index3j] - rhoinit[(index3i - 1) * NX + index3j]) * u[index3i * NX + index3j] +
			phi2 * (rhoinit[index3i * NX + index3j + 1] - rhoinit[index3i * NX + index3j - 1]) * v[index3i * NX + index3j]
			);
		initJd4 = abs(rhodiff[index4i * NX + index4j] + 
			phi1 * rhoinit[(index4i) * NX + index4j] * (u[(index4i + 1) * NX + index4j] - u[(index4i - 1) * NX + index4j]) +
			phi1 * rhoinit[index4i * NX + index4j] * (v[index4i * NX + index4j + 1] - v[index4i * NX + index4j - 1]) + 
			phi2 * (rhoinit[(index4i + 1) * NX + index4j] - rhoinit[(index4i - 1) * NX + index4j]) * u[index4i * NX + index4j] +
			phi2 * (rhoinit[index4i * NX + index4j + 1] - rhoinit[index4i * NX + index4j - 1]) * v[index4i * NX + index4j]
			);
		initJd0 = abs(rhodiff[i * NX + j] +
			phi1 * rhoinit[i * NX + j] * (u[(i + 1) * NX + j] - u[(i - 1) * NX + j]) +
			phi1 * rhoinit[i * NX + j] * (v[i * NX + j + 1] - v[i * NX + j - 1]) +
			phi2 * (rhoinit[(i + 1) * NX + j] - rhoinit[(i - 1) * NX + j]) * u[i * NX + j] +
			phi2 * (rhoinit[i * NX + j + 1] - rhoinit[i * NX + j - 1]) * v[i * NX + j]);
			
		initJd = initJd1 + initJd2 + initJd3 + initJd4 + initJd0;
		
		// IDK why, but doing this with the 9-stencil or whatever, gives terrible results, OR: the results are preferentially diagonal... idk why
		
		init1dux = u[(index1i + 1) * NX + index1j] - u[(index1i - 1) * NX + index1j];
		init1dvx = v[(index1i + 1) * NX + index1j] - v[(index1i - 1) * NX + index1j];
		init1dvy = v[index1i * NX + index1j + 1] - v[index1i * NX + index1j - 1];
		init1duy = u[index1i * NX + index1j + 1] - u[index1i * NX + index1j - 1];
		
		init2dux = u[(index2i + 1) * NX + index2j] - u[(index2i - 1) * NX + index2j];
		init2dvx = v[(index2i + 1) * NX + index2j] - v[(index2i - 1) * NX + index2j];
		init2dvy = v[index2i * NX + index2j + 1] - v[index2i * NX + index2j - 1];
		init2duy = u[index2i * NX + index2j + 1] - u[index2i * NX + index2j - 1];
		
		init3dux = u[(index3i + 1) * NX + index3j] - u[(index3i - 1) * NX + index3j];
		init3dvx = v[(index3i + 1) * NX + index3j] - v[(index3i - 1) * NX + index3j];
		init3dvy = v[index3i * NX + index3j + 1] - v[index3i * NX + index3j - 1];
		init3duy = u[index3i * NX + index3j + 1] - u[index3i * NX + index3j - 1];
		
		init4dux = u[(index4i + 1) * NX + index4j] - u[(index4i - 1) * NX + index4j];
		init4dvx = v[(index4i + 1) * NX + index4j] - v[(index4i - 1) * NX + index4j];
		init4dvy = v[index4i * NX + index4j + 1] - v[index4i * NX + index4j - 1];
		init4duy = u[index4i * NX + index4j + 1] - u[index4i * NX + index4j - 1];
		
		initJr =  init1dux * init1dux + init1duy * init1duy + init1dvx * init1dvx + init1dvy * init1dvy
			+ init2dux * init2dux + init2duy * init2duy + init2dvx * init2dvx + init2dvy * init2dvy
			+ init3dux * init3dux + init3duy * init3duy + init3dvx * init3dvx + init3dvy * init3dvy
			+ init4dux * init4dux + init4duy * init4duy + init4dvx * init4dvx + init4dvy * init4dvy;
		
	//	initJc = (u[i * NX + j] - uc[i * NX + j]) * (u[i * NX + j] - uc[i * NX + j]) 
	//		+ (v[i * NX + j] - vc[i * NX + j]) * (v[i * NX + j] - vc[i * NX + j]);
	
	//	initJc1 = gamma * exp((u[i * NX + j] - uc[i * NX + j]) * (u[i * NX + j] - uc[i * NX + j]));
	//	initJc2 = gamma * exp((v[i * NX + j] - vc[i * NX + j]) * (v[i * NX + j] - vc[i * NX + j]));
		initJc1 = gamma * (u[i * NX + j] - uc[i * NX + j]) * (u[i * NX + j] - uc[i * NX + j]);
		initJc2 = gamma * (v[i * NX + j] - vc[i * NX + j]) * (v[i * NX + j] - vc[i * NX + j]);
		initJc = initJc1 + initJc2;
		
		epsinit = initJd + alpha * initJr + initJc;
	//	epsinit = initJd + alpha * initJr + gamma * exp(initJc);
	//	epsinit = initJd + alpha * initJr;
		
		uinc = dist(gen);
		vinc = dist(gen);
		
		potJd1 = abs(rhodiff[index1i * NX + index1j] + 
			phi1 * rhoinit[(index1i) * NX + index1j] * ((u[(index1i + 1) * NX + index1j] + uinc) - u[(index1i - 1) * NX + index1j]) +
			phi1 * rhoinit[index1i * NX + index1j] * (v[index1i * NX + index1j + 1] - v[index1i * NX + index1j - 1]) +
			phi2 * (rhoinit[(index1i + 1) * NX + index1j] - rhoinit[(index1i - 1) * NX + index1j]) * u[index1i * NX + index1j] +
			phi2 * (rhoinit[index1i * NX + index1j + 1] - rhoinit[index1i * NX + index1j - 1]) * v[index1i * NX + index1j]
			);
		potJd2 = abs(rhodiff[index2i * NX + index2j] + 
			phi1 * rhoinit[(index2i) * NX + index2j] * (u[(index2i + 1) * NX + index2j] - (u[(index2i - 1) * NX + index2j] + uinc)) +
			phi1 * rhoinit[index2i * NX + index2j] * (v[index2i * NX + index2j + 1] - v[index2i * NX + index2j - 1]) +
			phi2 * (rhoinit[(index2i + 1) * NX + index2j] - rhoinit[(index2i - 1) * NX + index2j]) * u[index2i * NX + index2j] +
			phi2 * (rhoinit[index2i * NX + index2j + 1] - rhoinit[index2i * NX + index2j - 1]) * v[index2i * NX + index2j]
			);		
		potJd3 = abs(rhodiff[index3i * NX + index3j] + 
			phi1 * rhoinit[(index3i) * NX + index3j] * (u[(index3i + 1) * NX + index3j] - u[(index3i - 1) * NX + index3j]) +
			phi1 * rhoinit[index3i * NX + index3j] * (v[index3i * NX + index3j + 1] - (v[index3i * NX + index3j - 1] + vinc)) + 
			phi2 * (rhoinit[(index3i + 1) * NX + index3j] - rhoinit[(index3i - 1) * NX + index3j]) * u[index3i * NX + index3j] +
			phi2 * (rhoinit[index3i * NX + index3j + 1] - rhoinit[index3i * NX + index3j - 1]) * v[index3i * NX + index3j]
			);
		potJd4 = abs(rhodiff[index4i * NX + index4j] + 
			phi1 * rhoinit[(index4i) * NX + index4j] * (u[(index4i + 1) * NX + index4j] - u[(index4i - 1) * NX + index4j]) +
			phi1 * rhoinit[index4i * NX + index4j] * ((v[index4i * NX + index4j + 1] + vinc) - v[index4i * NX + index4j - 1]) + 
			phi2 * (rhoinit[(index4i + 1) * NX + index4j] - rhoinit[(index4i - 1) * NX + index4j]) * u[index4i * NX + index4j] +
			phi2 * (rhoinit[index4i * NX + index4j + 1] - rhoinit[index4i * NX + index4j - 1]) * v[index4i * NX + index4j]
			);
		potJd0 = abs(rhodiff[i * NX + j] +
			phi1 * rhoinit[i * NX + j] * (u[(i + 1) * NX + j] - u[(i - 1) * NX + j]) +
			phi1 * rhoinit[i * NX + j] * (v[i * NX + j + 1] - v[i * NX + j - 1]) + 
			phi2 * (rhoinit[(i + 1) * NX + j] - rhoinit[(i - 1) * NX + j]) * (u[i * NX + j] + uinc) +
			phi2 * (rhoinit[i * NX + j + 1] - rhoinit[i * NX + j - 1]) * (v[i * NX + j] + vinc));
			
		potJd = potJd1 + potJd2 + potJd3 + potJd4 + potJd0;
			
		pot1dux = u[(index1i + 1) * NX + index1j] - u[(index1i - 1) * NX + index1j] + uinc;
		pot1dvx = v[(index1i + 1) * NX + index1j] - v[(index1i - 1) * NX + index1j] + vinc;
		pot1dvy = v[index1i * NX + index1j + 1] - v[index1i * NX + index1j - 1];
		pot1duy = u[index1i * NX + index1j + 1] - u[index1i * NX + index1j - 1];
		
		pot2dux = u[(index2i + 1) * NX + index2j] - u[(index2i - 1) * NX + index2j] - uinc;
		pot2dvx = v[(index2i + 1) * NX + index2j] - v[(index2i - 1) * NX + index2j] - vinc;
		pot2dvy = v[index2i * NX + index2j + 1] - v[index2i * NX + index2j - 1];
		pot2duy = u[index2i * NX + index2j + 1] - u[index2i * NX + index2j - 1];
		
		pot3dux = u[(index3i + 1) * NX + index3j] - u[(index3i - 1) * NX + index3j];
		pot3dvx = v[(index3i + 1) * NX + index3j] - v[(index3i - 1) * NX + index3j];
		pot3dvy = v[index3i * NX + index3j + 1] - v[index3i * NX + index3j - 1] - vinc;
		pot3duy = u[index3i * NX + index3j + 1] - u[index3i * NX + index3j - 1] - uinc;
		
		pot4dux = u[(index4i + 1) * NX + index4j] - u[(index4i - 1) * NX + index4j];
		pot4dvx = v[(index4i + 1) * NX + index4j] - v[(index4i - 1) * NX + index4j];
		pot4dvy = v[index4i * NX + index4j + 1] - v[index4i * NX + index4j - 1] + vinc;
		pot4duy = u[index4i * NX + index4j + 1] - u[index4i * NX + index4j - 1] + uinc;
		
		potJr =   pot1dux * pot1dux + pot1duy * pot1duy + pot1dvx * pot1dvx + pot1dvy * pot1dvy
			+ pot2dux * pot2dux + pot2duy * pot2duy + pot2dvx * pot2dvx + pot2dvy * pot2dvy
			+ pot3dux * pot3dux + pot3duy * pot3duy + pot3dvx * pot3dvx + pot3dvy * pot3dvy
			+ pot4dux * pot4dux + pot4duy * pot4duy + pot4dvx * pot4dvx + pot4dvy * pot4dvy;
		
	//	potJc = (u[i * NX + j] - uc[i * NX + j] + uinc) * (u[i * NX + j] - uc[i * NX + j] + uinc) 
	//		+ (v[i * NX + j] - vc[i * NX + j] + vinc) * (v[i * NX + j] - vc[i * NX + j] + vinc);
		
	//	potJc1 = gamma * exp((u[i * NX + j] - uc[i * NX + j] + uinc) * (u[i * NX + j] - uc[i * NX + j] + uinc));
	//	potJc2 = gamma * exp((v[i * NX + j] - vc[i * NX + j] + vinc) * (v[i * NX + j] - vc[i * NX + j] + vinc));
		potJc1 = gamma * (u[i * NX + j] - uc[i * NX + j] + uinc) * (u[i * NX + j] - uc[i * NX + j] + uinc);
		potJc2 = gamma * (v[i * NX + j] - vc[i * NX + j] + vinc) * (v[i * NX + j] - vc[i * NX + j] + vinc);
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
			u[i * NX + j] = u[i * NX + j] + uinc;
			v[i * NX + j] = v[i * NX + j] + vinc;
			numchanges = numchanges + 1;
		}
		else{
			if (probofmove > distribution3(generator)) {
				u[i * NX + j] = u[i * NX + j] + uinc;
				v[i * NX + j] = v[i * NX + j] + vinc;
				numchanges = numchanges + 1;
			}
		} 
		
	
	//

	/*
	if ((s % 10000) == 0) {				//Checking sum of energy
		trueenergy = 0.0;
	//	for (int m = 2; m < NX - 2; m++) {		
	//		for (int n = 2; n < NX - 2; n++) {
		for (int m = 3; m < NX - 3; m++) {		
			for (int n = 3; n < NX - 3; n++) {
			
				initJd = abs(rhodiff[m * NX + n] + 
					phi1 * rhoinit[m * NX + n] * (u[(m + 1) * NX + n] - u[(m - 1) * NX + n]) +
					phi1 * rhoinit[m * NX + n] * (v[m * NX + n + 1] - v[m * NX + n - 1]) +
					phi2 * (rhoinit[(m + 1) * NX + n] - rhoinit[(m - 1) * NX + n]) * u[m * NX + n] +
					phi2 * (rhoinit[m * NX + n + 1] - rhoinit[m * NX + n - 1]) * v[m * NX + n]);
				init1dux = u[(m + 1) * NX + n] - u[(m - 1) * NX + n];
				init1dvx = v[(m + 1) * NX + n] - v[(m - 1) * NX + n];
				init1dvy = v[m * NX + n + 1] - v[m * NX + n - 1];
				init1duy = u[m * NX + n + 1] - u[m * NX + n - 1];
				initJr = init1dux * init1dux + init1duy * init1duy + init1dvx * init1dvx + init1dvy * init1dvy;
				
			//	initJc = (u[m * NX + m] - uc[m * NX + n]) * (u[m * NX + n] - uc[m * NX + n]) 
			//	+ (v[m * NX + n] - vc[m * NX + n]) * (v[m * NX + n] - vc[m * NX + n]);
		
			//	initJc1 = gamma * exp((u[m * NX + n] - uc[m * NX + n]) * (u[m * NX + n] - uc[m * NX + n));
			//	initJc2 = gamma * exp((v[m * NX + n] - vc[m * NX + n]) * (v[m * NX + n] - vc[m * NX + n]));
				initJc1 = gamma * (u[m * NX + n] - uc[m * NX + n]) * (u[m * NX + n] - uc[m * NX + n]);
				initJc2 = gamma * (v[m * NX + n] - vc[m * NX + n]) * (v[m * NX + n] - vc[m * NX + n]);
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
	double bindist = 2.0;
	
	// first we make an interpolation of the vector field
        for (int i = 0; i < 128; i++) {
        	for (int j = 0; j < 128; j++) {
        	
        		pointx = i * 2;
        		pointy = j * 2;
        		
        		Px = 0.0;
        		Py = 0.0;
        		
        		sumweight = 0.0;
        		weight = 1.0;
        		
        		IntDensity = 0.0;
        	
        		for (int n = 0; n < (nx * nx); n++){ 
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
        		
        		cout << pointx << "\t" << pointy << "\t" << Px << "\t" << Py << "\t" << IntDensity << endl;
        		
        	}
        }
	


}	// --- end of program ---



