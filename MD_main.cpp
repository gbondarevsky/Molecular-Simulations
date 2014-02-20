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
void velScale();
float totvelocsq();

int main(){
	genCoords();
	//printCoords();
	initveloc();
	printVel();
	velScale();
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

	srand((unsigned)time(0));

	for(int j=0; j<215; j++){
	r1 = number();//Creates random numbers for the normal distributions
	r2 = number();
	r3 = number();
	r4 = number();
	r5 = number(); 
	r6 = number();
	velocx[j] = 14.378*sqrt(T)*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2); // Assigns a random velocity in a normal distribution
	velocy[j] = 14.378*sqrt(T)*sqrt(-2.0*log(r3))*cos(2.0*M_PI*r4); //14.378 is the sqrt(k/m) 
	velocz[j] = 14.378*sqrt(T)*sqrt(-2.0*log(r5))*cos(2.0*M_PI*r6);
	}

	for( int k=0; k<N; k++){
		totalx = totalx + velocx[k];// sums all of the velocities
		totaly = totaly + velocy[k];
		totalz = totalz + velocz[k];
	}
	if (totalx != 0){
		float correctionx = totalx/N; // checks and corrections velocity to make total momentum 0
		for (int l=0; l<N; l++){
			velocx[l] = velocx[l] - correctionx;
		}
	}

	if (totaly != 0){
		float correctiony = totaly/N;
		for(int m=0; m<N; m++){
			velocy[m] = velocy[m] - correctiony;
		}
	}

	if (totalz != 0){
		float correctionz = totalz/N;
		for(int n=0; n<N; n++){
			velocz[n] = velocz[n] - correctionz;
		}
	}
	
}

float number(){
	rand(); rand(); rand(); // Magic
	float r = (float(rand()) / float(RAND_MAX));// Returns a random number between 0 and 1
	return r;
}

float kintemp(){
	float c = 7.48e-6;// Prefactor
	float kintemp = c*totvelocsq();
	return kintemp; 
}

float totvelocsq(){
	float totvelocsq;
	for(int i=0; i<215; i++){
		totvelocsq = totvelocsq + velocx[i]*velocx[i] + velocy[i]*velocy[i] + velocz[i]*velocz[i];
		}
	return totvelocsq;
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

void velScale(){
	float Tt = kintemp();
	float deltaT = abs(T - Tt);
	for ( ; deltaT > 2; ){
		for( int i=0; i<N; i++){
			velocx[i] = ((T+100)/(Tt+100))*velocx[i]; //Makes the scaling go slower 
			velocy[i] = ((T+100)/(Tt+100))*velocy[i];
			velocz[i] = ((T+100)/(Tt+100))*velocz[i];
			}
		Tt = kintemp();
		deltaT = abs( T -Tt);
	}
}