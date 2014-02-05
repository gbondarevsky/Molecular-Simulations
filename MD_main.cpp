#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>

using namespace std;

//Global Variables

//Function Prototypes
float genCoords();
void printarray(float arr[][3], int rows, int columns);


int main(){
        printarray(genCoords(),216,3);
        return 0;
}

float genCoords(){
        float coords [216][3]; //Rows are each particle. Columns are x, y, z, coords
        for (int i = 0; i <= 215; i++){
                coords[i][0] = float(i%18)*1.75; //Fill all x coords 
        }
        for (int i = 0; i <= 215; i++){
                coords[i][1] = float(i%16)*1.75; //Fill y coords
        }
        for (int i = 0; i <= 215; i++){
                coords[i][2] = float(i%6)*1.75; //Fill z coords
        }
        return coords[215][2];
}

void printarray(float arr[][3],int rows,int columns){
        for (int j = 0; j < columns; j++){
                for (int i = 0; i < rows; i++){
                        cout << arr[i][j] << "\t";
                }
                cout << "\n";
        }
}
