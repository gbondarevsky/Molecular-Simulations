#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

//Constants
const float r = 15; //Distance from one particle to another.  We have to go redo the y and z directions at some point
const float rh = r/2.0;
const int N = 216; //Number of particles
const int Nmax = N/3; //Maximum number of particles per plane
const int xmax = 18; //Number of particles with unique x values in a single plane [Ask Gary].
const float dt = 1; //Time step
const float dt2 = 2*dt; //2*Time step
const float dtsq = dt*dt; //Time step squared
const float eps = 119.8*0.001380649; // epsilon
const float sig = 3.405; // sigma


//Global Variables
float coords[N][3];
float velocx[N];
float velocy[N];
float velocz[N];
float T = 20;
float rx[N];
float ry[N];
float rz[N];
float ax[N];
float ay[N];
float az[N];
float rxold[N];
float ryold[N];
float rzold[N];
float sumvsq; //v^2
float rij[N][N];
float xij[N][N];
float yij[N][N];
float zij[N][N];
float lj[N][N];


//Function Prototypes
int genCoords();
void initveloc();
float number();
float kintemp();
void printCoords();
void printVel();
void simulation();

int main(){
	genCoords();
	printCoords();
	initveloc();
	printVel();
	simulation();

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
	velocx[j] = 14.378*100000*sqrt(T)*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2); // Assigns a random velocity in a normal distribution
	velocy[j] = 14.378*100000*sqrt(T)*sqrt(-2.0*log(r3))*cos(2.0*M_PI*r4); //14.378 is the sqrt(k/m) 
	velocz[j] = 14.378*100000*sqrt(T)*sqrt(-2.0*log(r5))*cos(2.0*M_PI*r6);
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

void simulation(){

        for (int i = 0; i < 3; i++){
		if(i == 0){
			for (int j = 0; j < N; j++){
                       		rx[j]=coords[j][0] ;
               		}	
		}
		if(i == 1){
			for (int j = 0; j < N; j++){
                        	ry[j]=coords[j][1] ;
                	}
		}
		if(i == 2){
			for (int j = 0; j < N; j++){
                        	rz[j]=coords[j][2] ;
                	}
		}
        }
	
	for (int j = 0; j < N; j++){
        	rxold[j] = rx[j] - velocx[j] * dt  ;
        	ryold[j] = ry[j] - velocy[j] * dt  ;
        	rzold[j] = rz[j] - velocz[j] * dt  ;
       	}

	float totalx;
	float totaly;
	float totalz;
	float rxnewI;
	float rynewI;
	float rznewI;
	float vxI;
	float vyI;
	float vzI;

	for(int t=1 ; t < 5; t++){

		//Distance Matrix

		for(int i=0; i<N; i++){

			for(int j=0; j<i; j++){

				rij[i][j] = sqrt((rx[i]-rx[j])(rx[i]-rx[j]) + (ry[i]-ry[j])(ry[i]-ry[j]) + (rz[i]-rz[j])(rz[i]-rz[j])  );

			}

		}

		//LJ Potential - Uses parameters for Argon

		for(int i=0; i<N; i++){

			for(int j=0; j<i; j++){
				
				LJ[i][j] = 4*eps*(pow(sig/rij[i][j],12) - pow(sig/rij[i][j],6));
		
			}

		}

		//Cartesian Distances

		for(int i=0; i<N; i++){

			for(int j=0; j<i; j++){

				xij = rx[i]-rx[j];
				yij = ry[i]-ry[j];
				zij = rz[i]-rz[j];

			}

		}

		//Minimum Image

		//Forces

		//Accelerations
	
		for( int i=0; i<N; i++){

			rxnewI = 2.0 * rx[i] - rxold[i] + dtsq * ax[i];
			rynewI = 2.0 * ry[i] - ryold[i] + dtsq * ay[i];
			rznewI = 2.0 * rz[i] - rzold[i] + dtsq * az[i];
			vxI = ( rxnewI - rxold[i] ) / dt2;
			vyI = ( rynewI - ryold[i] ) / dt2;
			vzI = ( rznewI - rzold[i] ) / dt2;
			sumvsq = sumvsq + vxI * vxI + vyI * vyI + vzI * vzI;	
			totalx = totalx + vxI;
			totaly = totaly + vyI;
			totalz = totalz + vzI;
			rxold[i] = rx[i];
			ryold[i] = ry[i];
			rzold[i] = rz[i];
			rx[i] = rxnewI;
			ry[i] = rynewI;
			rz[i] = rznewI;

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
	float totvelocsq;

	for(int i=0; i<215; i++){
		totvelocsq = totvelocsq + velocx[i]*velocx[i] + velocy[i]*velocy[i] + velocz[i]*velocz[i];
		}
	float kintemp = c*totvelocsq/10000000000;
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
