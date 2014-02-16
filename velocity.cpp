#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <cstdlib>
using namespace std;

void initveloc();// function prototypes
float number();

float velocx[216];// Global Variables
float velocy[216];
float velocz[216];
float T = 20;


int main() {
	initveloc();

	for( int i=0; i<215; i++){
		cout << velocx[i] << "\t" << velocy[i] << "\t" << velocz[i] << "\n";
	}
	

}

void initveloc(){
	
	float r1, r2, r3, r4, r5, r6;//Variable Declarations
	float totalx, totaly, totalz;

	srand(time(NULL));

	for(int j=0; j<215; j++){
	r1 = number();//Creates random numbers for the normal distributions
	r2 = number();
	r3 = number();
	r4 = number();
	r5 = number(); 
	r6 = number();
	velocx[j] = 14.378*sqrt(T)*sqrt(-2*log(r1))*cos(2*M_PI*r2); // Assigns a random velocity in a normal distribution
	velocy[j] = 14.378*sqrt(T)*sqrt(-2*log(r3))*cos(2*M_PI*r4);//14.378 is the sqrt(k/m) 
	velocz[j] = 14.378*sqrt(T)*sqrt(-2*log(r5))*cos(2*M_PI*r6);
	}

	for( int k=0; k<215; k++){
		totalx = totalx + velocx[k];// sums all of the velocities
		totaly = totaly + velocy[k];
		totalz = totalz + velocz[k];
	}
	if (totalx != 0){
		float correctionx = totalx/216; // checks and corrections velocity to make total momentum 0
		for (int l=0; l<215; l++){
			velocx[l] = velocx[l] - correctionx;
		}
	}

	if (totaly != 0){
		float correctiony = totaly/216;
		for(int m=0; m<215; m++){
			velocy[m] = velocy[m] - correctiony;
		}
	}

	if (totalz != 0){
		float correctionz = totalz/216;
		for(int n=0; n<215; n++){
			velocz[n] = velocz[n] - correctionz;
		}
	}
	
}
float number(){
	double r = ((double) rand() / (RAND_MAX));// Returns a random number between 0 and 1
	return r;
}
