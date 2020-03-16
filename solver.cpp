#include "C:\Users\robbe\Desktop\Code\visit_writer.h"
#include <math.h>
#include <cmath>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <algorithm>

#define XMIN -2
#define YMIN -2
#define ZMIN -2
#define XMAX 2
#define YMAX 2
#define ZMAX 2
#define CELLSIZE 0.1
#define GRIDSIZE 40 //Should be (XMAX - XMIN) / CELLSIZE
#define NUMCELLS 64000 //Should be GRIDSIZE cubed
#define GRAVITY -9.8 // gravity acceleration m/s^2
#define SOUND 1450.0 // speed of sound in m/s
#define RHO_0 10000 //reference density of water in kg/m^3
#define P_0 101325 //reference pressure of water in Pa
#define DIFF 0.0000001 //diffusion magnitude

#define ALPHA_FLUID -0.01e2  //viscosity for fluid
#define ALPHA_BOUNDARY 500e-2 //viscosity for boundary (should be high to prevent penetration)
#define BDENSFACTOR 1.5 //Density is increased in boundary particles

#define C1 3.0e-0  //stress tensor constants for granular material
#define C2 3e5
#define C3 1.2e1
#define PHI 1.23 //friction angle (radians)
#define KC 15000 //cohesion

//#define GAMMASTRETCH 0.99
//#define GAMMACOMPRESS 0.1
//#define ALPHA 1.0

#define DT 0.0004 //Time step size



float cutoff = 0.06; //Interaction neighbourhood radius

//////USEFUL FUNCTIONS//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Kernel function
float kernel(float r) {
	if (r >= 0 && r <= cutoff) {
		return 1. / 3.14159 / (pow(cutoff, 3))*(1 - 3. / 2. * pow((r / cutoff), 2) + 3. / 4. * pow((r / cutoff), 3));
	}
	else if (r > cutoff && r < (2 * cutoff)) {
		return 1. / 3.14159 / (pow(cutoff, 3)) * 1 / 4. * pow(2 - (r / cutoff), 3);
	}
	else {
		return 0;
	}
}

float kernel_test(float r) {
	if (r >= 0 && r <= cutoff) {
		return 1. / 3.14159 / (pow(cutoff, 4))*(1 - 3. * pow((r / cutoff), 1) + 9. / 4. * pow((r / cutoff), 2));
	}
	else if (r > cutoff && r < (2 * cutoff)) {
		return -1. / 3.14159 / (pow(cutoff, 4)) * 1 / 2. * pow(2 - (r / cutoff), 2);
	}
	else {
		return 0;
	}
}

float kernel_derivative(float r) {
	if (r < cutoff) {
		return -45.0 / 3.14159 / pow(cutoff, 6)*pow((cutoff - r), 2);
	}
	else {
		return 0;
	}

}

//Dot product
inline float dot_prod(float x1, float y1, float z1, float x2, float y2, float z2) {
	return x1*x2+y1*y2+z1*z2;
}

//Cross products
inline float cross_prod_x(float x1, float y1, float z1, float x2, float y2, float z2) {
	return y1*z2 - z1*y2;
}

inline float cross_prod_y(float x1, float y1, float z1, float x2, float y2, float z2) {
	return -x1*z2 + z1*x2;
}

inline float cross_prod_z(float x1, float y1, float z1, float x2, float y2, float z2) {
	return x1*y2 - y1*x2;
}

//DEFINE PARTICLES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class Particle {

public:

	//Constructor
	Particle() {
		xcoord = ycoord = zcoord = 0.;
		xvel = yvel = zvel = 0.;
		xacc = yacc = 0.;
		zacc = GRAVITY;
		flag = false;
		boundary = false;
	}

	~Particle(){}

	//Constructor
	Particle(float x, float y, float z) {
		xcoord = x;
		ycoord = y;
		zcoord = z;
		xvel = yvel = zvel = 0;
		xacc = yacc = 0;
		zacc = GRAVITY;
		flag = false;
		boundary = false;
	}

	//Constructor for boundary particles
	Particle(float x, float y, float z, bool b) {
		xcoord = x;
		ycoord = y;
		zcoord = z;
		xvel = yvel = zvel = 0;
		xacc = yacc = zacc = 0;
		boundary = b;
		flag = false;
	}

	//Constructor
	Particle(float x, float y, float z, float vx, float vy, float vz) {
		xcoord = x;
		ycoord = y;
		zcoord = z;
		xvel = vx;
		yvel = vy;
		zvel = vz;
		xacc = yacc = 0;
		zacc = GRAVITY;
		flag = false;
		boundary = false;
	}

	//Coordinates
	float xcoord;
	float ycoord;
	float zcoord;

	//Velocity
	float xvel;
	float yvel;
	float zvel;

	//Acceleration
	float xacc;
	float yacc;
	float zacc;

	//Index for tracking
	int index;

	//Physical parameters
	float mass = 1;
	float dens = RHO_0; //Default density water
	float press = 0; //Default pressure of 101325 Pa (kg/s^2/m) = 1 atm
	float delpressz=0;
	float delpressy=0;
	float delpressx=0;
	float diffusionx = 0;
	float diffusiony = 0;
	float diffusionz = 0;
	//float springx = 0;
	//float springy = 0;
	//float springz = 0;
	float sigma = 0;
	//float springlength = 0.02;
	//float kspring = 0;

	float newdens=RHO_0; //Default density water
	float newpress=P_0; //Default pressure of 101325 Pa (kg/s^2/m) = 1 atm
	float newdelpressz;
	float newdelpressy;
	float newdelpressx;
	float newsigma=0;
	//float newspringlength = 0.06;

	float vel_grad[3][3] = {0}; //velocity gradient field
	float strain_rate[3][3] = { 0 }; 
	float stress_rate[3][3] = { 0 };
	float strain_rate_squared[3][3] = { 0 };
	float stress_tensor[3][3] = { 0 };
	float stress_tensor_squared[3][3] = { 0 };
	float stress_accel[3] = {0};


	bool boundary;
	bool solid = false;
	bool flag;


	void set_dens(float x) {
		dens = (x + kernel(0)) / 23.0 *(1 + float(boundary)*BDENSFACTOR)+9250;
	}

	void set_delpress(float x, float y, float z) {
		delpressx = x;
		delpressy = y;
		delpressz = z;
	}

	void set_coord(float x, float y, float z) {
		xcoord = x;
		ycoord = y;
		zcoord = z;
	}

	void set_vel(float x, float y, float z) {
		xvel = x;
		yvel = y;
		zvel = z;
	}

	void set_flag() {
		flag = true;
	}


	//Calculate distance between particles
	inline float distance(Particle P) {
		return sqrt(pow(rab_x(P),2) + pow(rab_y(P), 2) + pow(rab_z(P), 2));
	}

	float rab_x(Particle P) {
		return xcoord - P.xcoord;
	}

	float rab_y(Particle P) {
		return ycoord - P.ycoord;
	}

	float rab_z(Particle P) {
		return zcoord - P.zcoord;
	}

	float vab_x(Particle P) {
		return xvel - P.xvel;
	}

	float vab_y(Particle P) {
		return yvel - P.yvel;
	}

	float vab_z(Particle P) {
		return zvel - P.zvel;	
	}

	//Calculate point density from particle pair
	float density(Particle P) {
		return kernel(distance(P));
	}

	float kernel_derivative_x(Particle P) {
		return kernel_derivative(distance(P))* rab_x(P) / distance(P);
	}

	float kernel_derivative_y(Particle P) {
		return kernel_derivative(distance(P))* rab_y(P) / distance(P);
	}

	float kernel_derivative_z(Particle P) {
		return kernel_derivative(distance(P))* rab_z(P) / distance(P);
	}

	//Calculate contribution to del of pressure due from one neighbour particle
	float del_pressure_mag(Particle P) {
		return (P.press / pow(P.dens, 2) + press / pow(dens, 2))*kernel_derivative(distance(P)); //in direction r1-r2
	}

	float del_pressure_x(Particle P) {
		return ((P.press / pow(P.dens, 2) + press / pow(dens, 2) + calculate_sigma(P))*kernel_derivative_x(P)); //x-component
	}

	float del_pressure_y(Particle P) {
		return ((P.press / pow(P.dens, 2) + press / pow(dens, 2) + calculate_sigma(P))*kernel_derivative_y(P)); //y-component
	}

	float del_pressure_z(Particle P) {
		return ((P.press / pow(P.dens, 2) + press / pow(dens, 2) + calculate_sigma(P))*kernel_derivative_z(P)); //z-component
	}

	//calculate pressure from current density
	void calculate_pressure(void) {
		press = 1000*pow(SOUND,0)*RHO_0/7.0*(pow(dens / RHO_0, 7) - 1);   //b*[(rho/rho0)^gamma -1]     gamma = 7, ref density = 1, b = speed of sound in medium squared at reference density
			//press = 101325 * (pow(dens / RHO_0, 7));// *pow(10, 30);
	}

	float calculate_sigma(Particle P) {
		float d = dot_prod(vab_x(P), vab_y(P), vab_z(P), rab_x(P), rab_y(P), rab_z(P));     //gamma*(mean speed of sound)cutoff*(v_1-v_2 dot r_1-r_2)/((r_1-r_2)^2+0.01 cutoff^2) / (mean density)
		float d2 = pow(distance(P),2);
		return (ALPHA_FLUID * SOUND * cutoff * (d / (d2 + 0.01*pow(cutoff, 2))) / ((dens + P.dens) / 2.0)) *(d < 0)*(1+!boundary*P.boundary*ALPHA_BOUNDARY);
		// alpha = 0.01 but may need to be tuned depending on cutoff

	}

	void update(void) {

		if (!flag) { //check if particle has already been updated

			set_dens(newdens);
			calculate_pressure();
			set_delpress(newdelpressx, newdelpressy, newdelpressz);

			for (int p = 0; p < 3; p++) {
				for (int q = 0; q < 3; q++) {
					stress_tensor[p][q] = DT*stress_rate[p][q];
				}
			}

			if (!boundary) {

				xcoord = xcoord + DT*xvel + DIFF*diffusionx;
				ycoord = ycoord + DT*yvel + DIFF*diffusiony;
				zcoord = zcoord + DT*zvel + DIFF*diffusionz;

				xvel = (xvel + DT*xacc + DT*(stress_accel[0]));
				yvel = (yvel + DT*yacc + DT*(stress_accel[1]));
				zvel = (zvel + DT*zacc + DT*(stress_accel[2]));

				//Acceleration due to grav
				xacc = -(300.0 / dens)*delpressx;
				yacc = -(300.0 / dens)*delpressy;
				zacc = GRAVITY + (-300.0 / dens)*delpressz;
			}
		}
		flag = false; //reset flag for next timestep
	}
};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


//DEFINE NODES///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


struct node
{
	Particle data;
	node *next;
};

class list
{
//private:
	//node *head, *tail;
public:
	node *head, *tail;
	list()
	{
		head = NULL;
		tail = NULL;
	}

	void createnode(Particle &value)
	{
		node *temp = new node;
		temp->data = value;
		temp->next = NULL;
		if (head == NULL)
		{
			head = temp;
			tail = temp;
			temp = NULL;
		}
		else
		{
			tail->next = temp;
			tail = temp;
		}
	}

	void insert_start(Particle &value)
	{
		node *temp = new node;
		temp->data = value;
		temp->next = head;
		head = temp;
	}

	void insert_position(int pos, Particle value)
	{
		node *pre = new node;
		node *cur = new node;
		node *temp = new node;
		cur = head;
		for (int i = 1; i<pos; i++)
		{
			pre = cur;
			cur = cur->next;
		}
		temp->data = value;
		pre->next = temp;
		temp->next = cur;
	}

	void delete_first()
	{
		node *temp = new node;
		temp = head;
		head = head->next;
		delete temp;
	}

	void delete_last()
	{
		node *current = new node;
		node *previous = new node;
		current = head;
		while (current->next != NULL)
		{
			previous = current;
			current = current->next;
		}
		tail = previous;
		previous->next = NULL;
		delete current;
	}

	void delete_position(int pos)
	{
		node *current = new node;
		node *previous = new node;
		current = head;
		
		if (pos == 1) {
			head = head->next;
			delete current;
		}

		for (int i = 1; i<pos; i++)
		{
			previous = current;
			current = current->next;
		}
		previous->next = current->next;
	}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Set up storage
const int npts = 1;  //number of fluid particles
const int nbpts = 1313;  //number of boundary particles
const int nspts = 1000; //number of solid particles

const int tpts = 10000;  //number of time points

const int ncoord = (npts) * 3; //number of coordinates
const int nbcoord = nbpts * 3; //number of coordinates
const int nscoord = nspts * 3; //number of coordinates

float sneighbour[nspts][nspts]; //keep track of distances between solid particles

int main(int argc, char **argv)
{

	float pts[ncoord]; //Array of particle coordinates for visit_writer
	float bpts[nbcoord]; //Array of boundary coordinates for visit_writer
	float spts[nscoord]; //Array of solid coordinates for visit_writer	

	list cells[NUMCELLS]; //set up a grid

	//Set up Solid Particles
	Particle *SPptr[nspts];
	for (int i = 0; i < nspts; i++) {
		for (int j = 0; j < nspts; j++) {
			if (i == 0) {
				SPptr[j] = new Particle(-.16 + 0.04*((j / 10) % 10), -0.16 + 0.04*(j / 10 / 10), -0.20+ (j % 10)*0.04, 0., 0., 0.);
				SPptr[j]->index = j;
				SPptr[j]->solid = true;
			} 
		}
	}

	//Set up Fluid Particles
	Particle *Pptr[npts+nbpts];
	for (int i = 0; i < npts; i++) {
		Pptr[i] = new Particle(0.04+0.04*((i/15)%15), 0.04*(i / 225), 0.0+(i%15)*0.04, 0., 0., 0.); 
	}
	
	//Set up boundary particles
	for (int i = 0; i < nbpts; i++) {
		Pptr[npts + i] = new Particle(-0.74 + 0.06*(i % 36), -0.74 + 0.06*(i / 36), -0.24, true);
	}

	//Set up initial nodes
	for (int i = 0; i < NUMCELLS; i++) {
		//fluid and boundary
		for (int j = 0; j < npts+nbpts; j++) {			
			if (Pptr[j]->xcoord > XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE && Pptr[j]->xcoord <= XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE + CELLSIZE) {          //x is int(i/(GRIDSIZE*GRIDSIZE)), y is int(i/GRIDSIZE), z is i%GRIDSIZE
				if (Pptr[j]->ycoord > YMIN + int((i%(GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE && Pptr[j]->ycoord <= YMIN + int((i % (GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE + CELLSIZE) {
					if (Pptr[j]->zcoord > ZMIN + (i%GRIDSIZE)*CELLSIZE && Pptr[j]->zcoord <= ZMIN + (i%GRIDSIZE)*CELLSIZE + CELLSIZE) {
						cells[i].createnode(*Pptr[j]);
						//std::cout << "Created a node. \n";
					}
				}
			}
		}
		//solid particles
		for (int k = 0; k < nspts; k++) {
			if (SPptr[k]->xcoord > XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE && SPptr[k]->xcoord <= XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE + CELLSIZE) {          //x is int(i/(GRIDSIZE*GRIDSIZE)), y is int(i/GRIDSIZE), z is i%GRIDSIZE
				if (SPptr[k]->ycoord > YMIN + int((i % (GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE && SPptr[k]->ycoord <= YMIN + int((i % (GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE + CELLSIZE) {
					if (SPptr[k]->zcoord > ZMIN + (i%GRIDSIZE)*CELLSIZE && SPptr[k]->zcoord <= ZMIN + (i%GRIDSIZE)*CELLSIZE + CELLSIZE) {
						cells[i].createnode(*SPptr[k]);
					}
				}
			}
		}
	}


	//Storage for output
	int vardims[] = { 1 , 1};   // Two scalars
	int vardims2[] = { 1 };   // One scalar
	int vardims3[] = { 1,1,1};
	
	float a[npts];
	float b[npts];
	float a2[nbpts];
	float a3[nspts];
	float b3[nspts];
	float c3[nspts];

	const char * const varnames[] = { "density", "zacc"};
	const char * const varnames2[] = { "density"};
	const char * const varnames3[] = { "stress_accel_x","stress_accel_y","stress_accel_z" };

	float *arrays[] = { (float*)a, (float*)b };
	float *arrays2[] = { (float*)a2};
	float *arrays3[] = { (float*)a3, (float*)b3, (float*)c3 };
	
	for (int t = 0; t < tpts; t++) {  //time steps
		std::cout << "t= " << t << "\n";
		int idx = 0;
		int idx2 = 0;
		int idx3 = 0;
		auto start = std::chrono::high_resolution_clock::now();

		/////////Search for neighbours
		for (int i = 0; i < NUMCELLS; i++) {
			node *N = cells[i].head; //check current cell
			//int pos = 1;
			while (N != NULL) {
				
				float tempdens = 0; 
				float tempdelpressx = 0;
				float tempdelpressy = 0;
				float tempdelpressz = 0;
				float tempsigma = 0;
				float tempdiffusionx = 0;
				float tempdiffusiony = 0;
				float tempdiffusionz = 0;
				
				N->node::data.vel_grad[0][0] = N->node::data.vel_grad[0][1] = N->node::data.vel_grad[0][2] = 0;
				N->node::data.vel_grad[1][0] = N->node::data.vel_grad[1][1] = N->node::data.vel_grad[1][2] = 0;
				N->node::data.vel_grad[2][0] = N->node::data.vel_grad[2][1] = N->node::data.vel_grad[2][2] = 0;
				N->node::data.stress_accel[0] = N->node::data.stress_accel[1] = N->node::data.stress_accel[2] = 0;

				for (int j = i-1; j <= i+1; j++) {
					for (int k = -GRIDSIZE; k <= GRIDSIZE; k += GRIDSIZE) {
						for (int l = -GRIDSIZE*GRIDSIZE; l <= GRIDSIZE*GRIDSIZE; l += GRIDSIZE*GRIDSIZE) {
							if (j + k + l >= 0 && j + k + l < NUMCELLS &&  i+j>= 0 && i + k >= 0 && i + l >= 0) {
								node *N2 = cells[j + k + l].head; //compare with all neighbour cells
								while (N2 != NULL) {
									float ds = N->node::data.distance(N2->node::data);
									if (ds <= (2*cutoff) && ds > 0) {

										float k = kernel(ds);
										float rabx = N->node::data.rab_x(N2->node::data);
										float raby = N->node::data.rab_y(N2->node::data);
										float rabz = N->node::data.rab_z(N2->node::data);
										float vabx = N->node::data.vab_x(N2->node::data);
										float vaby = N->node::data.vab_y(N2->node::data);
										float vabz = N->node::data.vab_z(N2->node::data);
										float dkx = kernel_derivative(ds)*rabx/ds;
										float dky = kernel_derivative(ds)*raby/ds;
										float dkz = kernel_derivative(ds)*rabz/ds;
										
										float dkxtest = kernel_test(ds)*rabx / ds;
										float dkytest = kernel_test(ds)*raby / ds;
										float dkztest = kernel_test(ds)*rabz / ds;

										float d = dot_prod(vabx, vaby, vabz, rabx, raby, rabz);
										float d2 = pow(ds, 2);
										float s = (ALPHA_FLUID * SOUND * (cutoff * (d / (d2 + 0.01*pow(cutoff, 2))) +50*1.0/SOUND*pow(cutoff * (d / (d2 + 0.01*pow(cutoff, 2))),2))/ ((N->node::data.dens + N2->node::data.dens) / 2.0)) *(d < 0)*(1 + (!N->node::data.boundary)*(N2->node::data.boundary) * ALPHA_BOUNDARY);
										
										float dpx = (N2->node::data.press / pow(N2->node::data.dens, 2) + N->node::data.press / pow(N->node::data.dens, 2) + s)*dkx;
										float dpy = (N2->node::data.press / pow(N2->node::data.dens, 2) + N->node::data.press / pow(N->node::data.dens, 2) + s)*dky;
										float dpz = (N2->node::data.press / pow(N2->node::data.dens, 2) + N->node::data.press / pow(N->node::data.dens, 2) + s)*dkz;
										
										//if (N->node::data.solid && N2->node::data.solid) {
											N->node::data.vel_grad[0][0] += -vabx*dkxtest / N2->node::data.dens;
											N->node::data.vel_grad[0][1] += -vaby*dkxtest / N2->node::data.dens;
											N->node::data.vel_grad[0][2] += -vabz*dkxtest / N2->node::data.dens;
											N->node::data.vel_grad[1][0] += -vabx*dkytest / N2->node::data.dens;
											N->node::data.vel_grad[1][1] += -vaby*dkytest / N2->node::data.dens;
											N->node::data.vel_grad[1][2] += -vabz*dkytest / N2->node::data.dens;
											N->node::data.vel_grad[2][0] += -vabx*dkztest / N2->node::data.dens;
											N->node::data.vel_grad[2][1] += -vaby*dkztest / N2->node::data.dens;
											N->node::data.vel_grad[2][2] += -vabz*dkztest / N2->node::data.dens;

											N->node::data.stress_accel[0] +=  (N->node::data.stress_tensor[0][0] * dkxtest + N->node::data.stress_tensor[0][1] * dkytest + N->node::data.stress_tensor[0][2] * dkztest) / pow(N->node::data.dens, 2) + (N2->node::data.stress_tensor[0][0] * dkxtest + N2->node::data.stress_tensor[0][1] * dkytest + N2->node::data.stress_tensor[0][2] * dkztest) / pow(N2->node::data.dens, 2);
											N->node::data.stress_accel[1] += (N->node::data.stress_tensor[1][0] * dkxtest + N->node::data.stress_tensor[1][1] * dkytest + N->node::data.stress_tensor[1][2] * dkztest) / pow(N->node::data.dens, 2) + (N2->node::data.stress_tensor[1][0] * dkxtest + N2->node::data.stress_tensor[1][1] * dkytest + N2->node::data.stress_tensor[1][2] * dkztest) / pow(N2->node::data.dens, 2);
											N->node::data.stress_accel[2] +=  (N->node::data.stress_tensor[2][0] * dkxtest + N->node::data.stress_tensor[2][1] * dkytest + N->node::data.stress_tensor[2][2] * dkztest) / pow(N->node::data.dens, 2) + (N2->node::data.stress_tensor[2][0] * dkxtest + N2->node::data.stress_tensor[2][1] * dkytest + N2->node::data.stress_tensor[2][2] * dkztest) / pow(N2->node::data.dens, 2);

									//	}
										tempdens += k*(1 + float(!N->node::data.boundary)*float(N2->node::data.boundary)*BDENSFACTOR);
										tempdelpressx += dpx;
										tempdelpressy += dpy;
										tempdelpressz += dpz;
										tempdiffusionx += 1 / N2->node::data.dens*dkx;
										tempdiffusiony += 1 / N2->node::data.dens*dky;
										tempdiffusionz += 1 / N2->node::data.dens*dkz;
									}
									N2 = N2->next;
								}
							}
						}
					}
				}
				////////////////Neighbour check done////////////////////
				///////////////Calculate final density and pressure///////////
				N->node::data.newdens = (tempdens);
				N->node::data.newdelpressx = tempdelpressx;
				N->node::data.newdelpressy = tempdelpressy;
				N->node::data.newdelpressz = tempdelpressz;
				N->node::data.newsigma = tempsigma;
				N->node::data.diffusionx = tempdiffusionx;
				N->node::data.diffusiony = tempdiffusiony;
				N->node::data.diffusionz = tempdiffusionz;
				
				if (N->node::data.solid) {
					float tr = 0; //trace of strain rate
					float tr2 = 0; //trace of stress tensor
					float tr3 = 0; //trace of stress tensor squared
					float tr4 = 0; //trace of stress tensor times strain rate
					float tr5 = 0; //trace of strain rate squared
					for (int p = 0; p < 3; p++) {
						for (int q = 0; q < 3; q++) {
							N->node::data.strain_rate[p][q] = 0.5*(N->node::data.vel_grad[p][q] + N->node::data.vel_grad[q][p]);
							N->node::data.stress_tensor_squared[p][q] = (N->node::data.stress_tensor[p][0] * N->node::data.stress_tensor[0][q] + N->node::data.stress_tensor[p][1] * N->node::data.stress_tensor[1][q] + N->node::data.stress_tensor[p][2] * N->node::data.stress_tensor[2][q]);
							N->node::data.strain_rate_squared[p][q] = (N->node::data.strain_rate[p][0] * N->node::data.strain_rate[0][q] + N->node::data.strain_rate[p][1] * N->node::data.strain_rate[1][q] + N->node::data.strain_rate[p][2] * N->node::data.strain_rate[2][q]);
							tr4 += N->node::data.stress_tensor[p][q] * N->node::data.strain_rate[q][p];
						}
						tr += N->node::data.strain_rate[p][p];
						tr2 += N->node::data.stress_tensor[p][p];
						tr3 += N->node::data.stress_tensor_squared[p][p];
						tr5+= N->node::data.strain_rate_squared[p][p];
						
					}
					float J2 = 0.5*(tr2 + tr3);


					for (int p = 0; p < 3; p++) {
						for (int q = 0; q < 3; q++) {
							N->node::data.stress_rate[p][q] = 3 * C1*N->node::data.press*(N->node::data.strain_rate[p][q] -1. / 3.*tr*(p == q))+ C1*C2*(tr4 + tr*N->node::data.press)/(pow(N->node::data.press,2) + 1e5)*N->node::data.stress_tensor[p][q] - C1*C3*(tr+tr5)*N->node::data.stress_tensor[p][q];
							if (3*tan(PHI) / (sqrt(9 + 12 * pow(tan(PHI),2)))*N->node::data.press + KC / (sqrt(9 + 12 * pow(tan(PHI), 2))) < J2 && J2!=0 && N->node::data.press>0) {
								N->node::data.stress_tensor[p][q] *= (3 * tan(PHI) / (sqrt(9 + 12 * pow(tan(PHI), 2)))*N->node::data.press + KC / (sqrt(9 + 12 * pow(tan(PHI), 2)))) / J2;
							}
						}
					}
				}
				N = N->next;
			}
		}
		
		////////////////////Update positions////////////////////////////
		float kin = 0;
		for (int i = 0; i < NUMCELLS; i++) {
			node *N = cells[i].head; //check current cell
			int pos = 1;
			while (N != NULL) {
				kin = std::max(kin,sqrt(pow(N->node::data.xvel, 2) + pow(N->node::data.yvel, 2) + pow(N->node::data.zvel, 2)));
				if (!N->node::data.flag && !N->node::data.boundary && !N->node::data.solid) {
					pts[(3 * idx)] = N->node::data.xcoord;
					pts[(3 * idx) + 1] = N->node::data.ycoord;
					pts[(3 * idx) + 2] = N->node::data.zcoord;
					a[idx] = (N->node::data).dens;
					b[idx] = (N->node::data).zacc;
					idx += 1;
				}
				if (!N->node::data.flag && N->node::data.boundary) {
					bpts[(3 * idx2)] = N->node::data.xcoord;
					bpts[(3 * idx2) + 1] = N->node::data.ycoord;
					bpts[(3 * idx2) + 2] = N->node::data.zcoord;
					a2[idx2] = (N->node::data).dens;
					//b2[idx2] = (N->node::data).zacc;
					idx2 += 1;
				}

				if (!N->node::data.flag && N->node::data.solid) {
					spts[(3 * idx3)] = N->node::data.xcoord;
					spts[(3 * idx3) + 1] = N->node::data.ycoord;
					spts[(3 * idx3) + 2] = N->node::data.zcoord;
					a3[idx3] = (N->node::data).stress_accel[0];
					b3[idx3] = (N->node::data).stress_accel[1];
					c3[idx3] = (N->node::data).stress_accel[2];
					idx3 += 1;
				}

				
				
				N->node::data.update();
				
				//Debugging
				/*if (N->node::data.index == nspts - 1) {
					std::cout << "Density = " << N->node::data.dens << "\n";
					std::cout << "Pressure = " << N->node::data.press << "\n";
					std::cout << "Velocity = (" << N->node::data.xvel << ", " << N->node::data.yvel << ", " << N->node::data.zvel << ")\n";
					std::cout << "Velocity_x gradient = (" << N->node::data.vel_grad[0][0] << ", " << N->node::data.vel_grad[1][0] << ", " << N->node::data.vel_grad[2][0] << ")\n";
					std::cout << "Velocity_y gradient = (" << N->node::data.vel_grad[0][1] << ", " << N->node::data.vel_grad[1][1] << ", " << N->node::data.vel_grad[2][1] << ")\n";
					std::cout << "Velocity_z gradient = (" << N->node::data.vel_grad[0][2] << ", " << N->node::data.vel_grad[1][2] << ", " << N->node::data.vel_grad[2][2] << ")\n";
					std::cout << "Strain rate = (" << N->node::data.strain_rate[0][0] << ", " << N->node::data.strain_rate[1][0] << ", " << N->node::data.strain_rate[2][0] << ")\n";
					std::cout << " (" << N->node::data.strain_rate[0][1] << ", " << N->node::data.strain_rate[1][1] << ", " << N->node::data.strain_rate[2][1] << ")\n";
					std::cout << " (" << N->node::data.strain_rate[0][2] << ", " << N->node::data.strain_rate[1][2] << ", " << N->node::data.strain_rate[2][2] << ")\n";
					std::cout << "Stress tensor = (" << N->node::data.stress_tensor[0][0] << ", " << N->node::data.stress_tensor[1][0] << ", " << N->node::data.stress_tensor[2][0] << ")\n";
					std::cout << " (" << N->node::data.stress_tensor[0][1] << ", " << N->node::data.stress_tensor[1][1] << ", " << N->node::data.stress_tensor[2][1] << ")\n";
					std::cout << " (" << N->node::data.stress_tensor[0][2] << ", " << N->node::data.stress_tensor[1][2] << ", " << N->node::data.stress_tensor[2][2] << ")\n";
					std::cout << "Stress Acceleration = (" << N->node::data.stress_accel[0] << ", " << N->node::data.stress_accel[1] << ", " << N->node::data.stress_accel[0] << ")\n";

				}*/
				
				
				////////////////Check if particles have moved to another cell/////////////////////
				float xlowb = XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE;
				float xupb = XMIN + int(i / (GRIDSIZE*GRIDSIZE))*CELLSIZE + CELLSIZE;
				float ylowb = YMIN + int((i % (GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE;
				float yupb = YMIN + int((i % (GRIDSIZE*GRIDSIZE)) / (GRIDSIZE))*CELLSIZE + CELLSIZE;
				float zlowb = ZMIN + (i%GRIDSIZE)*CELLSIZE;
				float zupb = ZMIN + (i%GRIDSIZE)*CELLSIZE + CELLSIZE;
				int newidx = i;

				if (N->node::data.xcoord < xlowb && N->node::data.xcoord > XMIN)
					newidx -= GRIDSIZE*GRIDSIZE;
				if (N->node::data.ycoord < ylowb && N->node::data.ycoord > YMIN)
					newidx -= GRIDSIZE;
				if (N->node::data.zcoord < zlowb && N->node::data.zcoord > ZMIN)
					newidx -= 1;
				if (N->node::data.xcoord > xupb && N->node::data.xcoord < XMAX)
					newidx += GRIDSIZE*GRIDSIZE;
				if (N->node::data.ycoord > yupb && N->node::data.ycoord < YMAX)
					newidx += GRIDSIZE;
				if (N->node::data.zcoord > zupb && N->node::data.zcoord < ZMAX)
					newidx += 1;
				
				if (newidx == i) {
					N = N->next;
					pos++;
				}
				else {
					cells[newidx].insert_start(N->node::data);
					if (newidx > i)
						cells[newidx].head->data.set_flag();
					N = N->next;
					cells[i].delete_position(pos);
				}
			}
		}
		
		if (t % 5 == 0) {
			//Write to file
			std::ostringstream oss;
			oss << "C:\\Users\\robbe\\Desktop\\Code\\anim_" << t/5 << ".vtk";
			std::string var = oss.str();
			const char* cstr = var.c_str();
			//write_point_mesh(cstr, 0, npts, pts, 2, vardims, varnames, arrays);
		}

		if (t % 10 == 0) {
			//Write each frame to file
			std::ostringstream oss;
			oss << "C:\\Users\\robbe\\Desktop\\Code\\anim_s" << t / 10 << ".vtk";
			std::string var = oss.str();
			const char* cstr = var.c_str();
			write_point_mesh(cstr, 0, nspts, spts, 3, vardims3, varnames3, arrays3);
		}

		if (t == 1) {
			//Write each frame to file
			std::ostringstream oss;
			oss << "C:\\Users\\robbe\\Desktop\\Code\\anim_boundary" << ".vtk";
			std::string var = oss.str();
			const char* cstr = var.c_str();
			write_point_mesh(cstr, 0, nbpts, bpts, 1, vardims2, varnames2, arrays2);
		}

		std::cout << "Max Kinetic Energy = " << kin << "\n\n";  //stability check
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = (stop - start)/1e9;
		std::cout << "Execution time: "<< duration.count() << "\n";
	}



	return 0;
}