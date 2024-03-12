#include "metropolis_hastings.h"
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

int seed1=21;
int nx=196;

// these are 1D, how to treat nx if nx != ny
std::default_random_engine generator;
	
std::normal_distribution<double> distribution(0.0,0.1);
std::uniform_real_distribution<double> distribution2(-1.0,1.0);
std::uniform_real_distribution<double> distribution3(0.0,1.0);
std::uniform_int_distribution<> distribution4(0, nx-1);
 	
std::uniform_int_distribution<> distributionx(2, nx-2);
std::uniform_int_distribution<> distributiony(2, nx-2);
 	
generator.seed(seed1);
 	
// delete the below and use: 
// std::normal_distribution<double> dist(0.0,0.05);
// If you can't install boost
 	
boost::mt19937 gen;
boost::random::normal_distribution<> dist(0.0,0.05);
gen.seed(seed1);
dist.reset();


void metropolis_hastings(double epsinit, double epsfinal, double epsdiff, int nosteps, double alpha, double temp, double phi1, double phi2, double gamma, double gamma2, double rhoinit[], double rhodiff[], double u[], double v[], double uc[], double vc[]){

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
		
		initJd1 = abs(rhodiff[index1i * nx + index1j] + 
			phi1 * rhoinit[(index1i) * nx + index1j] * (u[(index1i + 1) * nx + index1j] - u[(index1i - 1) * nx + index1j]) +
			phi1 * rhoinit[index1i * nx + index1j] * (v[index1i * nx + index1j + 1] - v[index1i * nx + index1j - 1]) +
			phi2 * (rhoinit[(index1i + 1) * nx + index1j] - rhoinit[(index1i - 1) * nx + index1j]) * u[index1i * nx + index1j] +
			phi2 * (rhoinit[index1i * nx + index1j + 1] - rhoinit[index1i * nx + index1j - 1]) * v[index1i * nx + index1j]
			);
		initJd2 = abs(rhodiff[index2i * nx + index2j] + 
			phi1 * rhoinit[(index2i) * nx + index2j] * (u[(index2i + 1) * nx + index2j] - u[(index2i - 1) * nx + index2j]) +
			phi1 * rhoinit[index2i * nx + index2j] * (v[index2i * nx + index2j + 1] - v[index2i * nx + index2j - 1]) +
			phi2 * (rhoinit[(index2i + 1) * nx + index2j] - rhoinit[(index2i - 1) * nx + index2j]) * u[index2i * nx + index2j] +
			phi2 * (rhoinit[index2i * nx + index2j + 1] - rhoinit[index2i * nx + index2j - 1]) * v[index2i * nx + index2j]
			);	
		initJd3 = abs(rhodiff[index3i * nx + index3j] + 
			phi1 * rhoinit[(index3i) * nx + index3j] * (u[(index3i + 1) * nx + index3j] - u[(index3i - 1) * nx + index3j]) +
			phi1 * rhoinit[index3i * nx + index3j] * (v[index3i * nx + index3j + 1] - v[index3i * nx + index3j - 1]) + 
			phi2 * (rhoinit[(index3i + 1) * nx + index3j] - rhoinit[(index3i - 1) * nx + index3j]) * u[index3i * nx + index3j] +
			phi2 * (rhoinit[index3i * nx + index3j + 1] - rhoinit[index3i * nx + index3j - 1]) * v[index3i * nx + index3j]
			);
		initJd4 = abs(rhodiff[index4i * nx + index4j] + 
			phi1 * rhoinit[(index4i) * nx + index4j] * (u[(index4i + 1) * nx + index4j] - u[(index4i - 1) * nx + index4j]) +
			phi1 * rhoinit[index4i * nx + index4j] * (v[index4i * nx + index4j + 1] - v[index4i * nx + index4j - 1]) + 
			phi2 * (rhoinit[(index4i + 1) * nx + index4j] - rhoinit[(index4i - 1) * nx + index4j]) * u[index4i * nx + index4j] +
			phi2 * (rhoinit[index4i * nx + index4j + 1] - rhoinit[index4i * nx + index4j - 1]) * v[index4i * nx + index4j]
			);
		initJd0 = abs(rhodiff[i * nx + j] +
			phi1 * rhoinit[i * nx + j] * (u[(i + 1) * nx + j] - u[(i - 1) * nx + j]) +
			phi1 * rhoinit[i * nx + j] * (v[i * nx + j + 1] - v[i * nx + j - 1]) +
			phi2 * (rhoinit[(i + 1) * nx + j] - rhoinit[(i - 1) * nx + j]) * u[i * nx + j] +
			phi2 * (rhoinit[i * nx + j + 1] - rhoinit[i * nx + j - 1]) * v[i * nx + j]);
			
		initJd = initJd1 + initJd2 + initJd3 + initJd4 + initJd0;
		
		// IDK why, but doing this with the 9-stencil or whatever, gives terrible results, OR: the results are preferentially diagonal... idk why
		
		init1dux = u[(index1i + 1) * nx + index1j] - u[(index1i - 1) * nx + index1j];
		init1dvx = v[(index1i + 1) * nx + index1j] - v[(index1i - 1) * nx + index1j];
		init1dvy = v[index1i * nx + index1j + 1] - v[index1i * nx + index1j - 1];
		init1duy = u[index1i * nx + index1j + 1] - u[index1i * nx + index1j - 1];
		
		init2dux = u[(index2i + 1) * nx + index2j] - u[(index2i - 1) * nx + index2j];
		init2dvx = v[(index2i + 1) * nx + index2j] - v[(index2i - 1) * nx + index2j];
		init2dvy = v[index2i * nx + index2j + 1] - v[index2i * nx + index2j - 1];
		init2duy = u[index2i * nx + index2j + 1] - u[index2i * nx + index2j - 1];
		
		init3dux = u[(index3i + 1) * nx + index3j] - u[(index3i - 1) * nx + index3j];
		init3dvx = v[(index3i + 1) * nx + index3j] - v[(index3i - 1) * nx + index3j];
		init3dvy = v[index3i * nx + index3j + 1] - v[index3i * nx + index3j - 1];
		init3duy = u[index3i * nx + index3j + 1] - u[index3i * nx + index3j - 1];
		
		init4dux = u[(index4i + 1) * nx + index4j] - u[(index4i - 1) * nx + index4j];
		init4dvx = v[(index4i + 1) * nx + index4j] - v[(index4i - 1) * nx + index4j];
		init4dvy = v[index4i * nx + index4j + 1] - v[index4i * nx + index4j - 1];
		init4duy = u[index4i * nx + index4j + 1] - u[index4i * nx + index4j - 1];
		
		initJr =  init1dux * init1dux + init1duy * init1duy + init1dvx * init1dvx + init1dvy * init1dvy
			+ init2dux * init2dux + init2duy * init2duy + init2dvx * init2dvx + init2dvy * init2dvy
			+ init3dux * init3dux + init3duy * init3duy + init3dvx * init3dvx + init3dvy * init3dvy
			+ init4dux * init4dux + init4duy * init4duy + init4dvx * init4dvx + init4dvy * init4dvy;
		
	//	initJc = (u[i * nx + j] - uc[i * nx + j]) * (u[i * nx + j] - uc[i * nx + j]) 
	//		+ (v[i * nx + j] - vc[i * nx + j]) * (v[i * nx + j] - vc[i * nx + j]);
	
	//	initJc1 = gamma * exp((u[i * nx + j] - uc[i * nx + j]) * (u[i * nx + j] - uc[i * nx + j]));
	//	initJc2 = gamma * exp((v[i * nx + j] - vc[i * nx + j]) * (v[i * nx + j] - vc[i * nx + j]));
		initJc1 = gamma * (u[i * nx + j] - uc[i * nx + j]) * (u[i * nx + j] - uc[i * nx + j]);
		initJc2 = gamma * (v[i * nx + j] - vc[i * nx + j]) * (v[i * nx + j] - vc[i * nx + j]);
		initJc = initJc1 + initJc2;
		
		epsinit = initJd + alpha * initJr + initJc;
	//	epsinit = initJd + alpha * initJr + gamma * exp(initJc);
	//	epsinit = initJd + alpha * initJr;
		
		uinc = dist(gen);
		vinc = dist(gen);
		
		potJd1 = abs(rhodiff[index1i * nx + index1j] + 
			phi1 * rhoinit[(index1i) * nx + index1j] * ((u[(index1i + 1) * nx + index1j] + uinc) - u[(index1i - 1) * nx + index1j]) +
			phi1 * rhoinit[index1i * nx + index1j] * (v[index1i * nx + index1j + 1] - v[index1i * nx + index1j - 1]) +
			phi2 * (rhoinit[(index1i + 1) * nx + index1j] - rhoinit[(index1i - 1) * nx + index1j]) * u[index1i * nx + index1j] +
			phi2 * (rhoinit[index1i * nx + index1j + 1] - rhoinit[index1i * nx + index1j - 1]) * v[index1i * nx + index1j]
			);
		potJd2 = abs(rhodiff[index2i * nx + index2j] + 
			phi1 * rhoinit[(index2i) * nx + index2j] * (u[(index2i + 1) * nx + index2j] - (u[(index2i - 1) * nx + index2j] + uinc)) +
			phi1 * rhoinit[index2i * nx + index2j] * (v[index2i * nx + index2j + 1] - v[index2i * nx + index2j - 1]) +
			phi2 * (rhoinit[(index2i + 1) * nx + index2j] - rhoinit[(index2i - 1) * nx + index2j]) * u[index2i * nx + index2j] +
			phi2 * (rhoinit[index2i * nx + index2j + 1] - rhoinit[index2i * nx + index2j - 1]) * v[index2i * nx + index2j]
			);		
		potJd3 = abs(rhodiff[index3i * nx + index3j] + 
			phi1 * rhoinit[(index3i) * nx + index3j] * (u[(index3i + 1) * nx + index3j] - u[(index3i - 1) * nx + index3j]) +
			phi1 * rhoinit[index3i * nx + index3j] * (v[index3i * nx + index3j + 1] - (v[index3i * nx + index3j - 1] + vinc)) + 
			phi2 * (rhoinit[(index3i + 1) * nx + index3j] - rhoinit[(index3i - 1) * nx + index3j]) * u[index3i * nx + index3j] +
			phi2 * (rhoinit[index3i * nx + index3j + 1] - rhoinit[index3i * nx + index3j - 1]) * v[index3i * nx + index3j]
			);
		potJd4 = abs(rhodiff[index4i * nx + index4j] + 
			phi1 * rhoinit[(index4i) * nx + index4j] * (u[(index4i + 1) * nx + index4j] - u[(index4i - 1) * nx + index4j]) +
			phi1 * rhoinit[index4i * nx + index4j] * ((v[index4i * nx + index4j + 1] + vinc) - v[index4i * nx + index4j - 1]) + 
			phi2 * (rhoinit[(index4i + 1) * nx + index4j] - rhoinit[(index4i - 1) * nx + index4j]) * u[index4i * nx + index4j] +
			phi2 * (rhoinit[index4i * nx + index4j + 1] - rhoinit[index4i * nx + index4j - 1]) * v[index4i * nx + index4j]
			);
		potJd0 = abs(rhodiff[i * nx + j] +
			phi1 * rhoinit[i * nx + j] * (u[(i + 1) * nx + j] - u[(i - 1) * nx + j]) +
			phi1 * rhoinit[i * nx + j] * (v[i * nx + j + 1] - v[i * nx + j - 1]) + 
			phi2 * (rhoinit[(i + 1) * nx + j] - rhoinit[(i - 1) * nx + j]) * (u[i * nx + j] + uinc) +
			phi2 * (rhoinit[i * nx + j + 1] - rhoinit[i * nx + j - 1]) * (v[i * nx + j] + vinc));
			
		potJd = potJd1 + potJd2 + potJd3 + potJd4 + potJd0;
			
		pot1dux = u[(index1i + 1) * nx + index1j] - u[(index1i - 1) * nx + index1j] + uinc;
		pot1dvx = v[(index1i + 1) * nx + index1j] - v[(index1i - 1) * nx + index1j] + vinc;
		pot1dvy = v[index1i * nx + index1j + 1] - v[index1i * nx + index1j - 1];
		pot1duy = u[index1i * nx + index1j + 1] - u[index1i * nx + index1j - 1];
		
		pot2dux = u[(index2i + 1) * nx + index2j] - u[(index2i - 1) * nx + index2j] - uinc;
		pot2dvx = v[(index2i + 1) * nx + index2j] - v[(index2i - 1) * nx + index2j] - vinc;
		pot2dvy = v[index2i * nx + index2j + 1] - v[index2i * nx + index2j - 1];
		pot2duy = u[index2i * nx + index2j + 1] - u[index2i * nx + index2j - 1];
		
		pot3dux = u[(index3i + 1) * nx + index3j] - u[(index3i - 1) * nx + index3j];
		pot3dvx = v[(index3i + 1) * nx + index3j] - v[(index3i - 1) * nx + index3j];
		pot3dvy = v[index3i * nx + index3j + 1] - v[index3i * nx + index3j - 1] - vinc;
		pot3duy = u[index3i * nx + index3j + 1] - u[index3i * nx + index3j - 1] - uinc;
		
		pot4dux = u[(index4i + 1) * nx + index4j] - u[(index4i - 1) * nx + index4j];
		pot4dvx = v[(index4i + 1) * nx + index4j] - v[(index4i - 1) * nx + index4j];
		pot4dvy = v[index4i * nx + index4j + 1] - v[index4i * nx + index4j - 1] + vinc;
		pot4duy = u[index4i * nx + index4j + 1] - u[index4i * nx + index4j - 1] + uinc;
		
		potJr =   pot1dux * pot1dux + pot1duy * pot1duy + pot1dvx * pot1dvx + pot1dvy * pot1dvy
			+ pot2dux * pot2dux + pot2duy * pot2duy + pot2dvx * pot2dvx + pot2dvy * pot2dvy
			+ pot3dux * pot3dux + pot3duy * pot3duy + pot3dvx * pot3dvx + pot3dvy * pot3dvy
			+ pot4dux * pot4dux + pot4duy * pot4duy + pot4dvx * pot4dvx + pot4dvy * pot4dvy;
		
	//	potJc = (u[i * nx + j] - uc[i * nx + j] + uinc) * (u[i * nx + j] - uc[i * nx + j] + uinc) 
	//		+ (v[i * nx + j] - vc[i * nx + j] + vinc) * (v[i * nx + j] - vc[i * nx + j] + vinc);
		
	//	potJc1 = gamma * exp((u[i * nx + j] - uc[i * nx + j] + uinc) * (u[i * nx + j] - uc[i * nx + j] + uinc));
	//	potJc2 = gamma * exp((v[i * nx + j] - vc[i * nx + j] + vinc) * (v[i * nx + j] - vc[i * nx + j] + vinc));
		potJc1 = gamma * (u[i * nx + j] - uc[i * nx + j] + uinc) * (u[i * nx + j] - uc[i * nx + j] + uinc);
		potJc2 = gamma * (v[i * nx + j] - vc[i * nx + j] + vinc) * (v[i * nx + j] - vc[i * nx + j] + vinc);
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
			u[i * nx + j] = u[i * nx + j] + uinc;
			v[i * nx + j] = v[i * nx + j] + vinc;
			numchanges = numchanges + 1;
		}
		else{
			if (probofmove > distribution3(generator)) {
				u[i * nx + j] = u[i * nx + j] + uinc;
				v[i * nx + j] = v[i * nx + j] + vinc;
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
}
