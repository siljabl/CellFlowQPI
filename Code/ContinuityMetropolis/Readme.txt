
The code currently runs like this:

	1.	use "TestRead.py" on DIC data with the output "dx1.txt" and "dy1.txt"
	2.	use "TestRead.py" on well data with the output "image1.txt" and "image2.txt"
			the fact that both of these are hard coded should be changed
	3.	compile the MetropolisVer8.cpp code with "g++ -std=c++11 -O3 -fopenmp -o output.x MetropolisVer8.cpp"
			note, you must have boost installed. Alternatively change to a random number generator of your choice
	4.	run with "./output.x 21 100000000 1.0 70.0 1.0 1.0 1.0 0.0286 > testout.txt"

Good choices for dx1 and dy1 are dx6-15.tif and dy6-15.tif, also Well1-1_resc_reg3.tif for image1 and Well1-1_resc_reg16.tif for image2
Currently there is something weird with the normalisations I have been given.
