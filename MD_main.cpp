#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

//Constants
const float r = 3.5; //Distance from one particle to another
const float rh = r/2.0;
const int N = 216; //Number of particles
const int Nmax = N/3; //Maximum number of particles per plane
const int xmax = 18; //Number of particles with unique x values in a single plane [Ask Gary].

//Global Variables
float coords[N][3];
float velocx[N];
float velocy[N];
float velocz[N];
float T = 20;


//Function Prototypes
int genCoords();
void initveloc();
float number();
float kintemp();
void printCoords();
void printVel();

int main(){
	genCoords();
	printCoords();
	initveloc();
	printVel();
    return 0;
}

int genCoords(){
    for (int i = 0; i < N; i++){
        coords[i][0] = (float(i%xmax)*rh) + (float(i/Nmax)*rh); //Fill all x coords
    }
    for (int i = 0; i < N; i++){
    	if (i%2 == 0){coords[i][1] = (float(2*((i%Nmax)/xmax))*r) + (float(i/Nmax)*rh);}
    	else {coords[i][1] = (float(2*((i%Nmax)/xmax)+1)*r) + (float(i/Nmax)*rh);}
    }
    for (int i = 0; i < N; i++){
		coords[i][2] = float(i/Nmax)*r;
    }
        return 0;
}

void printCoords(){
	cout << N << "\n";
	cout << "#" << "\n";
	for (int j = 0; j < N; j++){
				cout << "Ar ";
                for (int i = 0; i < 3; i++){
                        cout << coords[j][i] << " ";
                }
                cout << "\n";
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

float kintemp(){
	float c = 7.48e-6;// Prefactor
	float totvelocsq;

	for(int i=0; i<215; i++){
		totvelocsq = totvelocsq + velocx[i]*velocx[i] + velocy[i]*velocy[i] + velocz[i]*velocz[i];
		}
	float kintemp = c*totvelocsq;
	return kintemp; 
}

void printVel(){
cout << "Printing the components of velocity\n";
	for( int i=0; i<215; i++){
		cout << velocx[i] << " | " << velocy[i] << " | " << velocz[i] << "\n";
	}
	float t = kintemp();
	

	cout << "\n";
	cout <<"The kinetic temperature is " << t << "\n";
}