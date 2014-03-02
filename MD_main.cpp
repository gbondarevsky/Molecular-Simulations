#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

//Constants
const double kb = 1.38065e-33; //A^2 kg /fs^2 / K  Adam: Trust me leave it like this for now
const double r = 3.45; //Distance from one particle to another.  We have to go redo the y and z directions at some point
const double rh = r/2.0;
const int N = 216; //Number of particles
const int Nmax = N/3; //Maximum number of particles per plane
const int xmax = 18; //Number of particles with unique x values in a single plane [Ask Gary].
const double dt = 0.0001; //Time step in femptoseconds
const double dt2 = 2*dt; //2*Time step
const double dtsq = dt*dt; //Time step squared
const double eps = 119.8*kb; //*0.001380649; // epsilon
const double sig = 3.405; // sigma
const double mAr = 39.9/6.02e23/1000; //Mass of an Ar atom in kg
const double V = 10659.0;//Volume of the Box
const double boxl = 10.0/3.0 * pow(V,1.0/3.0);//In Angströms
const double boxx = 9.0/20.0 * boxl;
const double boxy = 8.0/20.0 * boxl;
const double boxz = 3.0/20.0 * boxl;
const double rcut = 2.5*sig; //Cutoff distance
const double T = 119.8;


//Global Variables
double coords[N][3];
double velocx[N];
double velocy[N];
double velocz[N];
double rx[N];
double ry[N];
double rz[N];
double ax[N];
double ay[N];
double az[N];
double rxold[N];
double ryold[N];
double rzold[N];
double sumvsq; //v^2
double rij[N][N];//pair distances in terms of r
double dxij[N][N];
double dyij[N][N];
double dzij[N][N];
double LJ[N][N];//LJ potential
double Fx[N][N];//Forces
double Fy[N][N];//Forces
double Fz[N][N];//Forces
bool neighborlist[N][N]; //Neighbor List
double totLJ; //Total potential energy
double rmsvdt;
double p;
double totalE;
double KE;

//Function Prototypes
int genCoords();
void initveloc();
double number();
double kintemp(double v);
void printCoords();
void printVel();
void velScale();
float totvelocsq();
void simulation();
void cartDist();
double minimage(double x1,double x2, double boxl);
void distMat();
void LJpot();
void Forces();
void Acceleration();
void neighbor();
void pressure();

int main(){
	genCoords();
	//printCoords();
	initveloc();
    velScale();
	//printVel();
	printf( "Time (ps)         TotalE        PE            KE            kintemp       pressure \n");
	simulation();
	cout << " Forces\n";
	for (int i=0; i<10; i++){
		cout << "\n";
		for (int j=0; j<10; j++){
			printf( " %11e ", Fx[i][j]);
		}
	}
	cout << "\n\n rij \n";
	for (int i=0; i<10; i++){
		cout << "\n";
		for (int j=0; j<10; j++){
			printf( " %11e ", rij[i][j]);
		}
	}
	cout << "\n\n dxij\n";
	for (int i=0; i<10; i++){
		cout << "\n";
		for (int j=0; j<10; j++){
			printf( " %11e ", dxij[i][j]);
		}
	}
	cout << "\n\n dyij\n";
	for (int i=0; i<10; i++){
		cout << "\n";
		for (int j=0; j<10; j++){
			printf( " %11e ", dyij[i][j]);
		}
	}
	cout << "\n\n dzij\n";
	for (int i=0; i<10; i++){
		cout << "\n";
		for (int j=0; j<10; j++){
			printf( " %11e ", dzij[i][j]);
		}
	}
	cout << "\n\n Accelerations\n";
	for( int i=0; i<10; i++){
		printf(" %11e %11e %11e \n", ax[i], ay[i], az[i]);
		}

	return 0;
}

int genCoords(){
    for (int i = 0; i < N; i++){
        coords[i][0] = (double(i%xmax)*rh) + (double(i/Nmax)*rh); //Fill all x coords
    }
    for (int i = 0; i < N; i++){
    	if (i%2 == 0){coords[i][1] = (double(2*((i%Nmax)/xmax))*r) + (double(i/Nmax)*rh);}
    	else {coords[i][1] = (double(2*((i%Nmax)/xmax)+1)*r) + (double(i/Nmax)*rh);}
    }
    for (int i = 0; i < N; i++){
		coords[i][2] = double(i/Nmax)*r;
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
	
	double r1, r2, r3, r4, r5, r6;//Variable Declarations
	double totalx = 0, totaly = 0, totalz = 0;
	double prefactor = sqrt(kb/mAr);
	srand((unsigned)time(0));

	for(int j=0; j<N; j++){
	r1 = number();//Creates random numbers for the normal distributions
	r2 = number();
	r3 = number();
	r4 = number();
	r5 = number(); 
	r6 = number();
	velocx[j] = prefactor*sqrt(T)*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2); // Assigns a random velocity in a normal distribution
	velocy[j] = prefactor*sqrt(T)*sqrt(-2.0*log(r3))*cos(2.0*M_PI*r4); //14.378 is the sqrt(k/m) 
	velocz[j] = prefactor*sqrt(T)*sqrt(-2.0*log(r5))*cos(2.0*M_PI*r6);
	}

	for( int k=0; k<N; k++){
		totalx = totalx + velocx[k];// sums all of the velocities
		totaly = totaly + velocy[k];
		totalz = totalz + velocz[k];
	}
	if ((totalx > 0.00001) || (totalx < -0.00001)){
		double correctionx = totalx/N; // checks and corrections velocity to make total momentum 0
		for (int l=0; l<N; l++){
			velocx[l] = velocx[l] - correctionx;
		}
	}

	if ((totaly > 0.00001) || (totaly < -0.00001)){
		double correctiony = totaly/N;
		for(int m=0; m<N; m++){
			velocy[m] = velocy[m] - correctiony;
		}
	}

	if ((totalz > 0.00001) || (totalz < -0.00001)){
		double correctionz = totalz/N;
		for(int n=0; n<N; n++){
			velocz[n] = velocz[n] - correctionz;
		}
	}
}

void simulation(){
	//This converts coords file to cartesian positions
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
	
	//This predicts old positions
	for (int j = 0; j < N; j++){
        rxold[j] = rx[j] - velocx[j] * dt  ;
        ryold[j] = ry[j] - velocy[j] * dt  ;
        rzold[j] = rz[j] - velocz[j] * dt  ;
    }

	//Local variables used in Verlet algorithm
	double totalx;
	double totaly;
	double totalz;
	double rxnewI;
	double rynewI;
	double rznewI;
	double vxI;
	double vyI;
	double vzI;

	neighbor();

	//loop over time
	for(int t=1 ; t < 20000; t++){
		totLJ = 0;
		sumvsq  = 0;
        totalx = 0;
        totaly = 0;
        totalz = 0;
        vxI = 0;
        vyI = 0;
        vzI = 0;
		rmsvdt = sqrt(sumvsq/N)*dt;
		for(int i=0; i < N; i++){
			for(int j=0; j<i; j++){
				if((rmsvdt > rij[i][j]-rcut) && (rij[i][j] > rcut)){
					neighbor();
					break;
				}
			}
		}
		//Distance Matrix
            distMat();
		//LJ Potential - Uses parameters for Argon
            LJpot();
         //Get Distance in Cartesians
            cartDist();
		//Forces
	          Forces();
		//Accelerations
			Acceleration();
		//Verlet Algorithm
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
	pressure();
	KE = 0.5*sumvsq*mAr;
     	totalE = totLJ +KE;

	if (t%500 == 0){
        	printf("%6lf    %10lf   %10lf    %10lf    %10lf    %10e\n",(t*dt/1000), (totalE/kb/T), (totLJ/kb/T), (KE/kb/T), kintemp(sumvsq), p);
	}

	}
}

double number(){
	rand(); rand(); rand(); // Magic
	double r = (double(rand()) / double(RAND_MAX));// Returns a random number between 0 and 1
	return r;
}

double kintemp(double v){
	double c = mAr/N/kb/3;// Prefactor
	double kintemp = c*v;
	return kintemp; 
}

float totvelocsq(){
	float totvelocsq = 0;
	for(int i=0; i<N; i++){
		totvelocsq = totvelocsq + velocx[i]*velocx[i] + velocy[i]*velocy[i] + velocz[i]*velocz[i];
		}
	return totvelocsq;
}

void printVel(){
cout << "Printing the components of velocity\n";
	for( int i=0; i<N; i++){
	//	cout << velocx[i] << " | " << velocy[i] << " | " << velocz[i] << "\n";
		printf( "%10f %10f %10f \n", velocx[i], velocy[i], velocz[i]) ;
	}
	double t = kintemp(totvelocsq());

	cout << "\n";
	cout <<"The kinetic temperature is " << t << "\n";
}

void cartDist(){
//Cartesian Distances
	for(int i=0; i<N; i++){
		for(int j=0; j<i; j++){
      		dxij[i][j] = abs(minimage( rx[i], rx[j], boxx));
      		dyij[i][j] = abs(minimage( ry[i], ry[j], boxy));
     		dzij[i][j] = abs(minimage( rz[i], rz[j], boxz));
   		}
	}
}                
//Min Image
double minimage(double x1, double x2, double boxl){
    double dist =  x1 - x2;
    dist = dist - floor( dist/boxl + 0.5)*boxl;
    return dist;
}


//Distance Matrix
void distMat(){
	for(int i=0; i<N; i++){
		for(int j=0; j<i; j++){
			rij[i][j] = sqrt(pow(minimage(rx[i],rx[j],boxx),2) + pow(minimage(ry[i],ry[j],boxy),2) + pow(minimage(rz[i],rz[j],boxz),2));
		}
	}
}

void neighbor(){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){ 
			if(rij[i][j] < rcut)
				neighborlist[i][j] = true;
			else
				neighborlist[i][j] = false;
		}
	}
}

//LJ Potential
void LJpot(){
	for(int i=0; i<N; i++){
		for(int j=0; j<i; j++){
            if(rij[i][j] < rcut){
                LJ[i][j] = 4.0*eps*(pow(sig/rij[i][j],12) - pow(sig/rij[i][j],6));
                totLJ = totLJ + LJ[i][j];
                //	cout << totLJ << " ";
            }
		}
	}
}

//Forces - The expression is completely obvious and not something that you should probably ask Chad
void Forces(){	
    for(int j=0; j<N; j++){
    	for(int i=0; i<N; i++){
    		if(rij[i][j] < rcut){
    			if(i > j){
				Fx[i][j] = -12.0*eps/(pow(2.0,1.0/6.0)*sig)*(pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),13) - pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),7))*dxij[i][j]/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j]));
				Fy[i][j] = -12.0*eps/(pow(2.0,1.0/6.0)*sig)*(pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),13) - pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),7))*dyij[i][j]/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j]));
				Fz[i][j] = -12.0*eps/(pow(2.0,1.0/6.0)*sig)*(pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),13) - pow(pow(2.0,1.0/6.0)*sig/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j])),7))*dzij[i][j]/(sqrt(dxij[i][j]*dxij[i][j]+dyij[i][j]*dyij[i][j]+dzij[i][j]*dzij[i][j]));
				
				}
				else if(i==j){
				Fx[i][j] = 0;
				Fy[i][j] = 0;
				Fz[i][j] = 0;
				}
				else {
				Fx[i][j] = -Fx[j][i];
				Fy[i][j] = -Fy[j][i];
				Fz[i][j] = -Fz[j][i];
				}
			}
       }
	}		
}

void Acceleration(){
    for (int i=0; i<N; i++) {
        ax[i] = 0;
        ay[i] = 0;
        az[i] = 0;
    }
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			ax[i] = ax[i] + Fx[i][j]/mAr;
			ay[i] = ay[i] + Fy[i][j]/mAr;
			az[i] = az[i] + Fz[i][j]/mAr;
		}
	}
}
void pressure(){
	double virial = 0;
	for(int i=0; i<N; i++){
		for(int j=0; j<i; j++){
		coords[i][0] = minimage(rx[i], rx[j],boxx);
		coords[i][1] = minimage(ry[i], ry[j],boxy);
		coords[i][2] = minimage(rz[i], rz[j],boxz);
		}
	}
	for(int i=0; i<N; i++){
		for( int j=0; j<i; j++){
			virial = virial + coords[i][0]*Fx[i][j] + coords[i][1]*Fy[i][j] + coords[i][2]*Fz[i][j];//Virial is the dot product of rij and fij
		}
	}
	virial = virial/3/V;// Since we are using an attractive potential the virial must always be negative
	if((totLJ < 0) && (virial < 0)){p = (N*kb*T/V) - virial;}
	else if((totLJ < 0) && (virial > 0)){p = (N*kb*T/V) + virial;}
    else if((totLJ > 0) && (virial < 0)){p = (N*kb*T/V) + virial;}
	else if((totLJ > 0) && (virial > 0)){p = (N*kb*T/V) - virial;}
    else {p = (N*kb*T/V);}
}

void velScale(){
	double Tt = kintemp(totvelocsq());
	double deltaT = abs(T - Tt);
	for ( ; deltaT > 0.1; ){
		for( int i=0; i<N; i++){
			velocx[i] = ((T+1000)/(Tt+1000))*velocx[i]; //Makes the scaling go slower
			velocy[i] = ((T+1000)/(Tt+1000))*velocy[i];
			velocz[i] = ((T+1000)/(Tt+1000))*velocz[i];
			}
		Tt = kintemp(totvelocsq());
        cout << Tt << "\n";
		deltaT = abs( T -Tt);
	}
}
