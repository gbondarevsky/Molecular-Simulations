#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>

using namespace std;

//Constants
const float r = 3.5; //Distance from one particle to another
const float rh = r/2.0;
const int N = 216; //Number of particles
const int Nmax = N/3; //Maximum number of particles per plane
const int xmax = 18; //Number of particles with unique x values in a single plane [Ask Gary].

//Global Variables
float coords[N][3];

//Function Prototypes
int genCoords();


int main(){
	genCoords();
	cout << N << "\n";
	cout << "#" << "\n";
	for (int j = 0; j < N; j++){
				cout << "Ar ";
                for (int i = 0; i < 3; i++){
                        cout << coords[j][i] << " ";
                }
                cout << "\n";
        }
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
