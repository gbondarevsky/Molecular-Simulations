#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>

using namespace std;

//Global Variables
float coords[216][3];

//Function Prototypes
int genCoords();


int main(){
	genCoords();
	cout << 216 << "\n";
	cout << "#" << "\n";
	for (int j = 0; j < 216; j++){
				cout << "Ar ";
                for (int i = 0; i < 3; i++){
                        cout << coords[j][i] << " ";
                }
                cout << "\n";
        }
        return 0;
}

int genCoords(){
    for (int i = 0; i <= 215; i++){
        coords[i][0] = (float(i%18)*1.75) + (float(i/72)*1.75); //Fill all x coords
    }
    for (int i = 0; i <= 215; i++){
    	if (i%2 == 0){coords[i][1] = (float(2*((i%72)/18))*3.5) + (float(i/72)*1.75);}
    	else {coords[i][1] = (float(2*((i%72)/18)+1)*3.5) + (float(i/72)*1.75);}
    }
    for (int i = 0; i <= 215; i++){
		coords[i][2] = float(i/72)*3.5;
    }
        return 0;
}
