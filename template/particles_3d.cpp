//
//  particles_3d.cpp
//
//
//  Created by Karandeep  Singh on 27/09/15.

// C++ program for particles diffusing in bulk and attaching to receptors diffusing on the membrane. Particles do 3d Brownian motion and receptors 2d Brownian motion. At some point of time, they attach to each other and once they both are bound, they are at rest. There has to be minimum number of receptors bound to a particle in order that the system can overcome bending energy of the membrane and have enough binding energy so as to wrap the particle. Which means cooperative effects of receptors are important. Once there are minimum number of receptors bound to a particle, particle is considered to be a bound. This minimum number of bound receptors needed is a free parameter, can be played with in order to compare with experiments. Experimentally, timescale is in minutes. so characteristic time should be in minutes.
//
//

#include<iostream>
#include<stdlib.h>
#include<time.h>
#include <algorithm>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <syslog.h>
#include <unistd.h>
#include "MersenneTwister.h"

int i,j,k,l,t,m,mmin=5,bp,br,mmax=30;  //mmin = minimum number of receptors bound to a particle , bp = bound particles, br = bound receptors
const int N = 10;                //number of particles
const int NR = 2000;             // number of receptors
const double ap = 0.5, rc = 1.0, rc2 = rc*rc;   //ap = particle radius, rc = critical distance between two particles for repulsive force calculations
//const double rhop = 0.005;
const double ar = 0.05, rhor = 0.01745;  //ar = receptor radius, rhor = receptor density on membrane, np = bulk viscosity, nr= membrane viscosity
const double Dp = 1.0, Dr = 1.0;  // Dp = particle diffusion coefficient calculated using Stoke-Einstein's formula =2 here, Dr = receptor diffusion coefficient = 2 here
long double time_steps = 1e6;   //total number of time steps
long double equil = 5e4;    //total number of time steps for equilibriation
double gamma_f = 1.0;   //friction in the system
const double dt = 1E-3;  //time step
double T = 1.0;       //temperature
//const double L = cbrt((4*3.142*ap*ap*ap*N)/(3*rhop));
const double LR = sqrt((3.142*ar*ar*NR)/(rhor));  //box length calculated using definition of density of receptors = total area of all receptors/box area
const double L = LR;   
const double H = 60.0;   //height of box
const double V = L*L*H;  //volume of box
const double rhop = (4.0*3.142*ap*ap*ap*N)/(3.0*V);   //volume fraction of particles = total volume of particles/volume of box, should be less than 0.7
double apatch = LR*LR*24.54*24.54*1e-6; // in microns    for calculating concentration of particles to match with experimental data, for 27 nm particles
double afactor = 137./apatch;      
double vfactor = V*24.54*24.54*24.54*1e-21; // in ml
double r2 = 0.0, xr[3], d[2];   
double h, hcut = 0.1684;   //hcut is the cutoff height when particles see receptors and bind with them
double d2 = 0.0;
//double d2cut = 4*ap*ap;
double d2cut = (ap+ar)*(ap+ar);   //d2cut is cutoff distance between particles and receptors when receptors can bind to particles
long int fSamp = 1e3;   //sampling frequency 
double realtime; //real time in minutes for 27 nm particles

using namespace std;

MTRand randi;

//tom:
vector<int> occN(N,0);
///

vector<double> receptor(3,0);  //defining arrays using vectors
vector< vector<double> >receptors(NR,receptor);  

vector<double> particle(3,0);
vector< vector<double> >particles(N,particle);

vector<double> force(3,0);
vector< vector<double> >forces(N,force);

vector<double> vmdpart(3,0);
vector< vector<double> >vmdparts(N,vmdpart);

vector<double> vmdrec(3,0);
vector< vector<double> >vmdrecs(NR,vmdrec);

vector<bool> occupiedreceptor(NR,false);
vector<bool> occupiedparticle(N,false);

void initPositions() //initializing particles on a lattice, probably from Daan-Frenkel or some other book
{
 	
    ofstream xii_out("initpos.txt");
    // Find M large enough to fit npart atoms on an fcc lattice
    int M = 1;
 
  
    // For 3-D, find M
        while (4*M*M*M<N) {
            M++;
        }
 
    // Prepare a FCC cubic cell
    double a2 = L/M;       // define unit size
    double a2z = H/M;				// define unit size
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};
 
    // Populate cubic/planar lattice
    int n = 0;
    for (int xx = 0; xx < M; xx++){
        for (int yy = 0; yy < M; yy++){
            for (int zz = hcut; zz < M; zz++){   //starting z coordinates of particles above hcut so that particles are not already in the adhesion zone from the beginning itself
                for (int k = 0; k < 4; k++){
                    if (n < N){
                        particles[n][0] = (xx + xCell[k])*a2;    //coordinates of particles
                        particles[n][1] = (yy + yCell[k])*a2;
                        particles[n][2] = (zz + zCell[k])*a2z;
          		vmdparts[n][0] = particles[n][0];     //coordinates for visualization purposes (vmd)
    			vmdparts[n][1] = particles[n][1];
    			vmdparts[n][2] = particles[n][2];
		//	cout << "xi= "<< particles[n][0] << "\t" << "yi= " << particles[n][1] << "\t" << "zi= " << particles[n][2] << "\n";
			xii_out << "xi= "<< particles[n][0] << "\t" << "yi= " << particles[n][1] << "\t" << "zi= " << particles[n][2] << "\n";   //writing particle coordinates in initpos.txt  
                        n++;
                    }
                }
            }
        }
    }
}



/*void initPositions()  //initializing particles randomly in the bulk
{

    for(int i=0; i<N; i++){
	particles[i][0] = randi.rand(1.0)*L;
	particles[i][1] = randi.rand(1.0)*L;
	particles[i][2] = hcut+randi.rand(1.0)*(L-hcut);
    vmdparts[i][0] = particles[i][0];
    vmdparts[i][1] = particles[i][1];
    vmdparts[i][2] = particles[i][2];
  //	particles[i][0] -= L*floor(particles[i][0]/(L));
 // 	particles[i][1] -= L*floor(particles[i][1]/(L));
//	if(particles[i][2]<=0.) particles[i][2] = -particles[i][2];
//	if(particles[i][2]>=L) particles[i][2] = L-(particles[i][2]-L);
//	cout << "xi= "<< particles[i][0] << "\t" << "yi= " << particles[i][1] << "\t" << "zi= " << particles[i][2] << "\n";
	//xii_out << "xi= "<< particles[i][0] << "\t" << "yi= " << particles[i][1] << "\t" << "zi= " << particles[i][2] << "\n";   
	 }
}*/

void initPositionsReceptor() //initializing receptors randomly on the bottom surface of the box
{

	for(int i=0; i<NR; i++){  //loop going through all receptors
	receptors[i][0] = randi.rand(1.0)*LR;  //x coordinates
	receptors[i][1] = randi.rand(1.0)*LR;  //y coordinates
        receptors[i][2] = 0.0;  //z coordinates of receptors are zero
        vmdrecs[i][0] = receptors[i][0]; //vmd coordinates 
        vmdrecs[i][1] = receptors[i][1];
        vmdrecs[i][2] = receptors[i][2];
 	  //  receptors[i][0]-=LR*floor(receptors[i][0]/(LR));
 	  //  receptors[i][1]-=LR*floor(receptors[i][1]/(LR));
 	  //  receptors[i][2]=0.0;;
	}
}

// Soft potential
void energy(void) //soft repulsive potential between particles, they can't overlap. there is no potential between receptors, they can overlap and doesn't matter since they are too small
{

    for(int i=0;i<(N-1);i++){ //loop going through 1st particle
        for(int j=(i+1);j<N;j++){  //through second particle
            xr[0]=particles[i][0]-particles[j][0];  //x coordinate between 2 particles
            xr[1]=particles[i][1]-particles[j][1];  //y coordinate "    "   "
            xr[2]=particles[i][2]-particles[j][2];  //z coordinate "   "   "
            xr[0]=xr[0]-L*nearbyint((xr[0]/L));  //periodic boundary conditions
            xr[1]=xr[1]-L*nearbyint((xr[1]/L));
            xr[2]=xr[2]-H*nearbyint((xr[2]/H));
            if (xr[0] == 0 && xr[1] == 0 && xr[2] == 0){    //check for particle overlaps
                cout << "xi= "<< particles[i][0] << "\t" << "yi= " << particles[i][1] << "\t" << "zi= " << particles[i][2] << "\n";
                cout << "xj= "<< particles[j][0] << "\t" << "yj= " << particles[j][1] << "\t" << "zi= " << particles[j][2] << "\n";
                cout << i << "\t" << j << "\n";
                exit(1);
            }
            for (int k=0; k<3; k++){  //for 3 dimensions
                r2 += xr[k]*xr[k];   //r^2 = x^2+y^2+z^2
            }
            if(r2<rc2){  //if r^2 < rc^2
                r2 = sqrt(r2);  
                double ff = ((2*ap)-r2)/(r2+0.01);  //soft repulsive force, taken from Marchetti's paper
                forces[i][0] += ff*xr[0];   //calculation of x,y,z forces on a pair of particles
                forces[i][1] += ff*xr[1];
                forces[i][2] += ff*xr[2];
                forces[j][0] -= ff*xr[0];
                forces[j][1] -= ff*xr[1];
                forces[j][2] -= ff*xr[2];
            }
        }
}
    }

/*void binding(void) //limited number of receptors per particle
{
    for(int i=0;i<N;i++){
        h = vmdparts[i][2];
        if(occN[i] <= mmax && h<hcut){
            for(int j=0;j<NR;j++){
                if(occN[i]>=mmax){
                   // cout << occN[i] <<endl;
                    break;
                }
                d[0] = vmdrecs[j][0]-vmdparts[i][0];
                d[1] = vmdrecs[j][1]-vmdparts[i][1];
                d[0] -= LR*nearbyint(d[0]/LR);
                d[1] -= LR*nearbyint(d[1]/LR);
                d2 = (d[0]*d[0])+(d[1]*d[1]);
		if(d2<d2cut && occupiedreceptor[j]==false){
                    occupiedreceptor[j]=true;
                    occupiedparticle[i]=true;
                    occN[i]++;
//		cout << i << "\t" << j << "\t" << occupiedparticle[i] << "\t" << occupiedreceptor[j] << "\t" << occN[i] << endl;
                }
            }
	        bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);  //counting bound particles
    		br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true); //counting bound receptors	
		cout << i  << "\t" << bp << "\t" << br << "\t" << occN[i] << endl;
        }
   	  
     }
}*/

/*
void binding(void) //hard binding of particles to receptors. It means that there has to be a minimum number of receptors bound to a particle so that particle can be counted as bound.
{
    for(int i=0;i<N;i++){ //loop through particles
        h = vmdparts[i][2];  //z coordinates of particles
        cout << h << endl;
	if(h<hcut){   //if h<hcut
//            cout << h << endl;
	  for(int j=0;j<NR;j++){  //loop through all receptors
               // if(occN[i]>=mmax){
                   // cout << occN[i] <<endl;
                //    break;
               // }
	    d[0] = vmdrecs[j][0]-vmdparts[i][0];  //distance between receptors and particles
	    d[1] = vmdrecs[j][1]-vmdparts[i][1];
	    d[0] -= LR*nearbyint(d[0]/LR);  //periodic bc
	    d[1] -= LR*nearbyint(d[1]/LR);
	    d2 = (d[0]*d[0])+(d[1]*d[1]);  //d^2 = dx^2+dy^2
              // cout << d2 << endl;
	    if(d2<d2cut && occupiedreceptor[j]==false){ 
	      occupiedreceptor[j]==true;  //receptor is bound
	      occN[i]++;	  //occupation number for particles
	    }

	    if(occN[i]>=mmin && occupiedparticle[i]==false) occupiedparticle[i]==true; //if occupation number > minimum no. of receptors bound to a particle, particle is bound.
	    else occupiedreceptor[j]==false;  //receptor is unbound
	  }
	}
	//cout << occN[i] << endl;
    }
}*/

/*void binding(void) //hard binding of particles to receptors. It means that there has to be a minimum number of receptors bound to a particle so that particle can be counted as bound.
{
    for(int i=0;i<N;i++){ //loop through particles
        h = vmdparts[i][2];  //z coordinates of particles
      //  cout << h << endl;
	if(h<hcut){   //if h<hcut
         // cout << h << endl;
	  for(int j=0;j<NR;j++){  //loop through all receptors
	    d[0] = vmdrecs[j][0]-vmdparts[i][0];  //distance between receptors and particles
	    d[1] = vmdrecs[j][1]-vmdparts[i][1];
	    d[0] -= LR*nearbyint(d[0]/LR);  //periodic bc
	    d[1] -= LR*nearbyint(d[1]/LR);
	    d2 = (d[0]*d[0])+(d[1]*d[1]);  //d^2 = dx^2+dy^2
	    if(d2<d2cut && occupiedreceptor[j]==false && occN[i]<=mmax){ //if cutoff distance satisfied, free receptors, and occ number is less than maximum 
//	    if(d2<d2cut && occupiedreceptor[j]==false){
	      occN[i]++;	  //occupation number for particles
	//    cout << i << "\t" << j << "\t" << occN[i] << endl; 
	     } 
	  }
	  if(occN[i] >= mmin){  //if occupation number for particles is greater than minimum
	//    cout << i << "\t" << occN[i] << endl;
	    occupiedparticle[i]=true;  //particle is bound
	//   cout << i << "\t" << std::boolalpha << occupiedparticle[i] << endl;
	    for(int j=0;j<NR;j++){  //loop through all receptors
	      d[0] = vmdrecs[j][0]-vmdparts[i][0];  //distance between receptors and particles
	      d[1] = vmdrecs[j][1]-vmdparts[i][1];
	      d[0] -= LR*nearbyint(d[0]/LR);  //periodic bc
	      d[1] -= LR*nearbyint(d[1]/LR);
	      d2 = (d[0]*d[0])+(d[1]*d[1]);  //d^2 = dx^2+dy^2
	      if(d2<d2cut && occupiedreceptor[j]==false && occN[i]<=mmax){ //if receptors satisfy cutoff distance and they are unbound
//	      if(d2<d2cut && occupiedreceptor[j]==false){	
		occupiedreceptor[j]=true;  //receptor is bound
//		cout << i << "\t" << j << "\t" << occupiedparticle[i] << "\t" << occupiedreceptor[j] << endl;
	      }
	    }
	  }
	//	bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);  //counting bound particles
    	//	br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true); //counting bound receptors
	//	cout << i  << "\t" << bp << "\t" << br << "\t" << occN[i] << endl;

	  else {
	    occN[i] = 0;  //else particle is free and receptor is free
	  //  occupiedreceptor[j]=false;
	  //  occupiedparticle[i]=false;
		}
	//	bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);  //counting bound particles
    	//	br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true); //counting bound receptors
	//	cout << i  << "\t" << bp << "\t" << br << "\t" << occN[i] << endl;

	}
	//cout << occN[i] << endl;
    }
}*/

void binding(void) // unlimited number of receptors per particle
    {
        for(int i=0;i<N;i++){
            h = vmdparts[i][2];
            if(h<hcut){
//               cout << h << endl;
		 for(int j=0;j<NR;j++){
                    d[0] = vmdrecs[j][0]-vmdparts[i][0];
                    d[1] = vmdrecs[j][1]-vmdparts[i][1];
                    d[0] = d[0]-LR*nearbyint(d[0]/LR);
                    d[1] = d[1]-LR*nearbyint(d[1]/LR);
                    d2 = (d[0]*d[0])+(d[1]*d[1]);
                    if(d2<d2cut && occupiedreceptor[j]==false){
                        occupiedreceptor[j]=true;
                        occupiedparticle[i]=true;
                    }
                }
//		bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);  //counting bound particles
  //  		br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true); //counting bound receptors
//		cout << i  << "\t" << bp << "\t" << br << "\t" << occN[i] << endl;

            }
        }
}
    

	/*
        if(h<hcut){
           for(int j=0;j<NR;j++){
	     //tom:
            
               if (m>=mmax){
		 cout << m <<endl; 
		 break;
	       }
               d[0] = vmdrecs[j][0]-vmdparts[i][0];
               d[1] = vmdrecs[j][1]-vmdparts[i][1];
               d[0] -= LR*nearbyint(d[0]/LR);
               d[1] -= LR*nearbyint(d[1]/LR);
               d2 = (d[0]*d[0])+(d[1]*d[1]);
               if(d2<d2cut && occupiedreceptor[j]==false){
                  occupiedreceptor[j]=true;
                  occupiedparticle[i]=true;
                  m++;
               }
               //m = count(occupiedreceptor.begin(), occupiedreceptor.end(), true);
               //if(m>=mmax) break;
           }
           //if(m>=mmax) continue;
        }
        //if(m>=mmax) continue;
	*/

int main()  //main program starts
{
    //cout << L << "\n";

    unsigned int saat = (unsigned int)time(0);
    randi.seed(saat);  //random seed
    
    openlog("particles_3d.cpp", LOG_PID|LOG_CONS, LOG_USER);

    initPositions();  //initialize particles
    initPositionsReceptor();  //initialize receptors

    double sqdt = sqrt(dt);  
    for(t=1; t<=equil;t++){  //equilibriation loop for particles and receptors
        energy();  //repulsive forces acting between particles
        for(int i=0;i<N;i++){  //loop through particles
            particles[i][0] = particles[i][0] + ((forces[i][0]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt; //x coordinates. overdamped diffusion. simulating Langevin equation or Brownian dynamics.
            vmdparts[i][0] = particles[i][0] - L*floor((particles[i][0]/L));  //periodic bc, these coordinates are used in visualizing in vmd
            particles[i][1] = particles[i][1] + ((forces[i][1]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt;  //y coordinates
            particles[i][2] = particles[i][2] + ((forces[i][2]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt;  //z coordinates
            vmdparts[i][1] = particles[i][1] - L*floor((particles[i][1]/L));  //periodic bc
            if(particles[i][2]<=hcut) particles[i][2] = hcut-(particles[i][2]-hcut);  //reflecting bc in z direction. particles cannot enter at the bootom or come out of the bottom, because there is a membrane.
            if(particles[i][2]>=H) particles[i][2] = H-(particles[i][2]-H);
            vmdparts[i][2] = particles[i][2];  
        }

        for(int i=0;i<NR;i++){
            receptors[i][0] = receptors[i][0] +sqrt(2*Dr)*randi.randNorm(0,1)*sqdt; //x coordinates for receptors. there is no force term here.
            vmdrecs[i][0] = receptors[i][0] - LR*floor((receptors[i][0]/LR));   //periodic bc
            receptors[i][1] = receptors[i][1] +sqrt(2*Dr)*randi.randNorm(0,1)*sqdt;  //y coordinates
            vmdrecs[i][1] = receptors[i][1] - LR*floor((receptors[i][1]/LR));  //periodic bc
            receptors[i][2] = 0.0;  //z coordinates = 0
            vmdrecs[i][2] = 0.0;
        }

    }

    ofstream x_out("coordinates3d_soft.xyz");  //writing particle coordinates for vmd
 
    ofstream xi_out("coordinates3d_soft.txt");  //writing coordinates for analysis

    ofstream xR_out("coordinatesR2d_soft.xyz");   //writing receptor coordinates for vmd
  
    ofstream xRi_out("coordinatesR2d_soft.txt");  //writing receptor coordinates for analysis
   
   // ofstream xii_out("initpos.txt");

//    ofstream occupancy_out("occupiedparticles.txt", ios::app);

    ofstream coord_x_out("xdata3d.txt");  //particle coordinates
    ofstream coord_y_out("ydata3d.txt");
    ofstream coord_z_out("zdata3d.txt");

    ofstream coord_xR_out("xRdata2d.txt"); //receptor coordinates
    ofstream coord_yR_out("yRdata2d.txt");

    //test
    x_out << 10 << "\n";  //for vmd purposes
    x_out << "i\n";
    xR_out << 2000 << "\n"; //for vmd purposes
    xR_out << "\n";
    
    for(vector< vector<double> > ::iterator it=vmdparts.begin(); it!=vmdparts.end(); it++){ //writing particle coordinates for vmd
        x_out << "C\t" << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
        // cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
    }
    syslog(LOG_INFO, "frame written");
    
    for(vector< vector<double> > ::iterator it=vmdrecs.begin(); it!=vmdrecs.end(); it++){   //writing receptor coordinate for vmd
        xR_out << "N\t" << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
        // cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
    }
    
    
    for(vector< vector<double> > ::iterator it=particles.begin(); it!=particles.end(); it++){  //particle coordinates for analysis
        xi_out << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
        //   cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
    }
    
    for(vector< vector<double> > ::iterator it=receptors.begin(); it!=receptors.end(); it++){  //receptor coordinates for analysis
        xRi_out << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
        //  cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
    }
    
    for(vector< vector<double> > ::iterator it=particles.begin(); it!=particles.end(); it++){
        coord_x_out << it[0][0] << "\t";
        coord_y_out << it[0][1] << "\t";
        coord_z_out << it[0][2] << "\t";
    }
    for(vector< vector<double> > ::iterator it=receptors.begin(); it!=receptors.end(); it++){
        coord_xR_out << it[0][0] << "\t";
        coord_yR_out << it[0][1] << "\t";
    }
    
    coord_x_out << "\n";
    coord_y_out << "\n";
    coord_z_out << "\n";
    coord_xR_out << "\n";
    coord_yR_out << "\n";
    // cout << "time: " << t << endl;
    
    ofstream occupancy_out("occupiedparticles.txt", ios::app); //writing ouput file
    bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);  //counting bound particles
    br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true); //counting bound receptors
    realtime = (1e-6*13.5*13.5*t*dt)/(0.55*0.55); //real time in seconds  //D=l^2/t, texp/tsim = lexp^2*Dsim/lsim^2Dexp, lexp=13.5nm Dexp=1um^2/s lsim=0.55(ap+ar) Dsim=1, tsim=dt*t  
   // realtime = (93e-6*ap*Dp*t*dt)/60.; //real time in minutes for 27 nm particles
    occupancy_out << "#1:lengthofbox 2.height 3.no. of particles 4.real particle concentration 5.no. of receptors 6.particle volume fraction 7.receptor density 8.part diffusion const 9.rec diff const 10.min number of rec required 11.max no. of rec bound to a particle 12.sim time 13.real time 14.bound particles 15.bound receptors 16.bound part/rbc 17.bound rec/rbc" << "\n";
    occupancy_out << L << "\t" << H << "\t" << N << "\t" << N/vfactor << "\t" << NR << "\t" << rhop << "\t" << rhor << "\t" << Dp << "\t" << Dr << "\t" << mmin << "\t" << mmax << "\t" << 0*dt << "\t" << realtime << "\t" << bp << "\t" << br << "\t" << bp*afactor << "\t" << br*afactor << "\n";
    
//occupancy_out << "#1:lengthofbox 2.height 3.no. of particles 4.real particle concentration 5.no. of receptors 6.particle volume fraction 7.receptor density 8.part diffusion const 9.rec diff const 10.min number of rec required 11.max no. of rec bound to a particle 12.sim time 13.real time 14.bound particles 15.bound receptors 16.bound part/rbc 17.bound rec/rbc" << "\n";  //writing lots of things in output file for analysis. boxlenght, boxheight, number of particles, particle concentration in bulk, number of receptors, volume fraction of particles, density of receptors, particle diffusion coefficient, receptor diffusion coefficient, minimum bound receptors, time steps, bound particles, bound receptors, bound particles per cell, bound receptors per cell
    occupancy_out.close();
    
    //test
    
    for(t=1; t<=time_steps;t++){ //dynamics loop for particles and receptors
        energy(); //forces
        for(int i=0;i<N;i++){ //particles
            if(occupiedparticle[i]==false){ //this loop is only valid for unbound particles because bound particles are already at rest. rest everything is same as above.
              particles[i][0] = particles[i][0] + ((forces[i][0]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt;
              vmdparts[i][0] = particles[i][0] - L*floor((particles[i][0]/L));
              particles[i][1] = particles[i][1] + ((forces[i][1]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt;
              vmdparts[i][1] = particles[i][1] - L*floor((particles[i][1]/L));
	      particles[i][2] = particles[i][2] + ((forces[i][2]*dt)/gamma_f)+sqrt(2*Dp)*randi.randNorm(0,1)*sqdt;
	      if(particles[i][2]<=0.) particles[i][2] = -particles[i][2];
              if(particles[i][2]>=H) particles[i][2] = H-(particles[i][2]-H);            	        
              vmdparts[i][2] = particles[i][2];
//	      cout << vmdparts[i][2] << endl;	
	    }
    	}
        
	for(int i=0;i<NR;i++){  //receptors
            if(occupiedreceptor[i]==false){ //only valid for unbound receptors.
            receptors[i][2] = 0.0;
            receptors[i][0] = receptors[i][0] +sqrt(2*Dr)*randi.randNorm(0,1)*sqdt;
            vmdrecs[i][0] = receptors[i][0] - LR*floor((receptors[i][0]/LR));
            receptors[i][1] = receptors[i][1] +sqrt(2*Dr)*randi.randNorm(0,1)*sqdt;
            vmdrecs[i][1] = receptors[i][1] - LR*floor((receptors[i][1]/LR));
            vmdrecs[i][2] = 0.0;
            }
        }
        binding();  //check for binding events

	if (!(t%fSamp)) {  //write data in files but after every fsamp steps

            x_out << 10 << "\n";
            x_out << "i\n";
            xR_out << 2000 << "\n";
            xR_out << "\n";

            for(vector< vector<double> > ::iterator it=vmdparts.begin(); it!=vmdparts.end(); it++){
                x_out << "C\t" << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
               // cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
            }
            syslog(LOG_INFO, "frame written");

            for(vector< vector<double> > ::iterator it=vmdrecs.begin(); it!=vmdrecs.end(); it++){
                xR_out << "N\t" << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
                // cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
            }


            for(vector< vector<double> > ::iterator it=particles.begin(); it!=particles.end(); it++){
                xi_out << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
             //   cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
            }

            for(vector< vector<double> > ::iterator it=receptors.begin(); it!=receptors.end(); it++){
                    xRi_out << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
                  //  cout << it[0][0] << "\t" <<  it[0][1] << "\t" << it[0][2] << "\n";
            }

            for(vector< vector<double> > ::iterator it=particles.begin(); it!=particles.end(); it++){
                coord_x_out << it[0][0] << "\t";
                coord_y_out << it[0][1] << "\t";
                coord_z_out << it[0][2] << "\t";
            }
            for(vector< vector<double> > ::iterator it=receptors.begin(); it!=receptors.end(); it++){
                coord_xR_out << it[0][0] << "\t";
                coord_yR_out << it[0][1] << "\t";
            }
           
            coord_x_out << "\n";
            coord_y_out << "\n";
            coord_z_out << "\n";
            coord_xR_out << "\n";
            coord_yR_out << "\n";
           // cout << "time: " << t << endl; 	
            
            ofstream occupancy_out("occupiedparticles.txt", ios::app);
            bp = count(occupiedparticle.begin(), occupiedparticle.end(), true);
    	    br = count(occupiedreceptor.begin(), occupiedreceptor.end(), true);
 	    realtime = (1e-6*13.5*13.5*t*dt)/(0.55*0.55);

           // realtime = (93e-6*ap*Dp*t*dt)/60.; //real time in minutes for 27 nm particles
            occupancy_out << L << "\t" << H << "\t" << N << "\t" << N/vfactor  << "\t" << NR << "\t" << rhop << "\t" << rhor << "\t" << Dp << "\t" << Dr << "\t" << mmin << "\t" << mmax << "\t" << t*dt << "\t" << realtime << "\t" << bp << "\t" << br << "\t" << bp*afactor << "\t" << br*afactor << "\n";
//	    occupancy_out << L << "\t" << H << "\t" << N << "\t"<< N/vfactor << "\t" << NR << "\t" << rhop << "\t" << rhor << "\t" << Dp << "\t" << Dr << "\t" << mmin << "\t" << mmax << "\t" << t*dt << "\t" << realtime << "\t" << bp << "\t" << br << "\t" << bp*afactor << "\n";
            occupancy_out.close();
        }
    }
    
    closelog ();
}
