#include "LObstacle.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// this function initializes all positions of the matrix to NULL
void init_matrix(Particle* matrix[][matrix_width]) {
	for (int i=0;i<matrix_height;i++)
		for (int j=0;j<matrix_width;j++)
			matrix[i][j]=NULL;
}

// this function initializes Totals
void init_totals(double totals[][8]) {
	for (int i=0;i<5;i++)
		for (int j=0;j<8;j++)
			totals[i][j]=0;
}

// This function checks the boundaries on the matrix
int check_boundaries(int final_axis_pos,int matrix_dim_size) {
	if (final_axis_pos<0)
		final_axis_pos=matrix_dim_size-abs(final_axis_pos);
	else if (final_axis_pos>matrix_dim_size-1)
		final_axis_pos=final_axis_pos-matrix_dim_size;
	return final_axis_pos;
}

// This function checks the matrix for erroneous obstacles
int check_matrix(int size,Particle* matrix[][matrix_width]) {
	Particle* part=NULL;
	LObstacle* lobst=NULL;
	int id=0,init_x=0,init_y=0,j=0,k=0,ecount=0;
	int posx,posy;
	bool res=1;
	
	for (int xxx=0;xxx<matrix_height;xxx++)
		for (int yyy=0;yyy<matrix_width;yyy++) {
			res=1;
			part=matrix[xxx][yyy];
			if ((part!=NULL) && (part->get_type()==4)) {
				lobst=(LObstacle*)part;
				id=lobst->get_id();
				init_x=lobst->get_init_obst_pos_x();
				init_y=lobst->get_init_obst_pos_y();
				j=0;
				while (res && (j<sqrt(size))) {
					k=0;
					while (res && (k<sqrt(size))) {
						posx=check_boundaries(init_x+j,matrix_height);
						posy=check_boundaries(init_y+k,matrix_width);
						if ((matrix[posx][posy]==NULL) || (matrix[posx][posy]->get_id()!=id))
							res=0;
						k++;
					}
					j++;
				}
			}
			if (res==0)
				ecount++;
		}
	return ecount;
}

// this function checks other position's availability
bool check_availability(int ri_h,int ri_w,int size,Particle* matrix[][matrix_width]) {
	bool res=1;
	int j=0,k=0;
	int posx,posy;
	
	while (res && (j<sqrt(size))) {
		k=0;
		while (res && (k<sqrt(size))) {
			posx=check_boundaries(ri_h+j,matrix_height);
			posy=check_boundaries(ri_w+k,matrix_width);
			if (matrix[posx][posy]!=NULL)
				res=0;
			k++;
		}
		j++;
	}
	return res;
}

// this function inserts particles at random positions in the matrix
void init_matrix_particles(int nparticles,int type,int size,double wfi,Particle* matrix[][matrix_width]) {
	int ri_h,ri_w;
	int posx,posy;
	bool create=1;
	
	int i=0;
	while (i<nparticles) {
		ri_h=rand()%(matrix_height); // produces a random number from 0 to matrix height -1
		ri_w=rand()%(matrix_width); // produces a random number from 0 to matrix_width -1
		create=check_availability(ri_h,ri_w,size,matrix); // checks if the matrix already contains a particle at the random position
		if (create) {
			for (int j=0;j<sqrt(size);j++) {
				for (int k=0;k<sqrt(size);k++) {
					posx=check_boundaries(ri_h+j,matrix_height);
					posy=check_boundaries(ri_w+k,matrix_width);
					if ((type==4) && (size>1)) {
						LObstacle* particle=new LObstacle(i+1,type,posx,posy,0,ri_h,ri_w);
						matrix[posx][posy]=particle;
					}
					else {
						Particle* particle=new Particle(i+1,type,posx,posy,0);
						particle->init_spin_particle();
						matrix[posx][posy]=particle;
						if ((type==4) && (((double)rand()/(RAND_MAX))<wfi))
							particle->set_interact(1);
						else 
							particle->set_interact(0);
					}
				}
			}
			i++;
		}
	}
}

// this function adds info about the particle
void add_particle_info(Particle* particle,int type,int iteration,double totals[][8]) {
	double lifetime=iteration-particle->get_it_creation();
	double moves=particle->get_part_move();
	double disp=particle->get_sq_disp();
	
	totals[0][type]+=1; // # of particles
	totals[1][type]=((totals[0][type]-1)*totals[1][type]+lifetime)/(totals[0][type]); // Lifetime
	totals[2][type]=((totals[0][type]-1)*totals[2][type]+moves)/(totals[0][type]); // Moves
	totals[3][type]=((totals[0][type]-1)*totals[3][type]+disp)/(totals[0][type]); // Square Displacement
	if (lifetime!=0) totals[4][type]=((totals[0][type]-1)*totals[4][type]+(disp/lifetime))/(totals[0][type]); // Diffusion
		
}

// this function writes the position of particles to a file
void write_particles_pos(int iteration,Particle* matrix[][matrix_width],FILE* file3) {
	Particle* particle;
	for (int i=0;i<matrix_height;i++) {
		for (int j=0;j<matrix_width;j++) {
			particle=matrix[i][j];
			if ((particle!=NULL) && (particle->get_type()!=4))
				fprintf(file3,"%i\t %i\t %i\t %i\t \n",iteration,particle->get_type(),i,j); // Writes spatial reaction position to file
		}
	}
}

// this function checks the number of particles in the matrix
void check_particles(Particle* matrix[][matrix_width],int* xxx) {
	Particle* particle;
	for (int i=0;i<matrix_height;i++)
		for (int j=0;j<matrix_width;j++) {
			particle=matrix[i][j];
			if (particle!=NULL)
				xxx[particle->get_type()]++;
		}
}

// this function converts char* to double
double to_double(const char *p) {
	std::stringstream ss(p);
	double result=0;
	ss>>result;
	return result;
}

// this is the main function
int main (int argc, char* const argv[]) {
	// Example on call - ./Diffusion Results1.txt Results2.txt Results3.txt Results4.txt 1 1 1 1000 0.01 0.1 0 4 0 1 0.02 0.04 0.01 0.01
	FILE* file1;
	FILE* file2;
	//FILE* file3;
	//FILE* file4;
	double spin,wfi,iterations,conc_A,conc_B,conc_O,size_O,mobil_O,f_rate,r_rate,c_rate,wf_rate,wr_rate;
	bool user_input=1;
	
	if (user_input) {
		file1=fopen(argv[1],"wb");
		file2=fopen(argv[2],"wb");
		//file3=fopen(argv[3],"wb");
		//file4=fopen(argv[4],"wb");
		spin=to_double(argv[5])/100;
		wfi=to_double(argv[6])/100;
		iterations=to_double(argv[7]);
		conc_A=to_double(argv[8])/100;
		conc_B=to_double(argv[9])/100;
		conc_O=to_double(argv[10])/100;
		size_O=to_double(argv[11]);
		mobil_O=to_double(argv[12]);
		f_rate=to_double(argv[13])/100;
		r_rate=to_double(argv[14])/100;
		c_rate=to_double(argv[15])/100;
		wf_rate=to_double(argv[16])/100;
		wr_rate=to_double(argv[17])/100;
	}
	else {
		file1=fopen("Results1.txt","wb");
		file2=fopen("Results2.txt","wb");
		//file3=fopen("Results3.txt","wb");
		//file4=fopen("Results4.txt","wb");
		spin=0;
		wfi=0;
		iterations=1000;
		conc_A=0.01;
		conc_B=0.1;
		conc_O=0.30;
		size_O=100;
		mobil_O=1;
		f_rate=1;
		r_rate=0.02;
		c_rate=0.04;
		wf_rate=1;
		wr_rate=1;
	}
	
	int matrix_size=matrix_height*matrix_width;
	double nparticles_A=floor(conc_A*matrix_size);
	double nparticles_B=floor(conc_B*matrix_size);
	double nparticles_O=floor(conc_O*matrix_size/size_O);
	if (size_O>1) wfi=0; // Weak force interactions may only occur if obstacles have the same size as any other particle
	Particle* matrix[matrix_height][matrix_width];
		
	// initializes the random seed and the matrix of particles
	srand((unsigned)time(NULL));
	init_matrix(matrix);
	init_matrix_particles(nparticles_O,4,size_O,wfi,matrix); // The identifier for the obstacle particles is 4!
	init_matrix_particles(nparticles_A,0,1,wfi,matrix);
	init_matrix_particles(nparticles_B,1,1,wfi,matrix);
	
	double totals[5][8]={0};
	int react[5]={0};
	int res[8]={0};
	int np_changes[8]={0};
	double nparticles[8]={0};
	double conc_particles[8]={0};
	double f_reactions=0,gamma=0,dist=0;
	Particle* particle;
	LObstacle* lobst;
	
	nparticles[0]=nparticles_A;
	nparticles[1]=nparticles_B;
	nparticles[4]=nparticles_O;
	conc_particles[0]=nparticles[0]/matrix_size;
	conc_particles[1]=nparticles[1]/matrix_size;
	conc_particles[4]=nparticles[4]*size_O/matrix_size;
	int total_nparticles=nparticles[0]+nparticles[1];
	if (mobil_O) total_nparticles+=nparticles[4];
	
	// MAIN CYCLE
	int i=1;
	while ((i<=iterations) && (conc_particles[1]>=0.01)) {
		init_totals(totals);
		f_reactions=0;
		for (int k=0;k<8;k++) np_changes[k]=0;
		
		int m=0;
		while (m<total_nparticles) {
			dist=0;
			int ri_h=rand()%(matrix_height); // produces a random number from 0 to matrix height -1
			int ri_w=rand()%(matrix_width); // produces a random number from 0 to matrix_width -1
			particle=matrix[ri_h][ri_w];
			if ((particle!=NULL) && (particle->get_it_move()!=i)) {
				if ((particle->get_type()==4) && (size_O>1) && mobil_O) {
					lobst=(LObstacle*)particle;
					if ((lobst->get_pos_x()==(lobst->get_init_obst_pos_x())) && (lobst->get_pos_y()==lobst->get_init_obst_pos_y())) {
						lobst->move(size_O,matrix);
						add_particle_info(lobst,lobst->get_type(),i,totals);
						m++;
					}
				}
				else {
					if ((particle->get_type()!=4) || mobil_O) {
						particle->move(i,spin,mobil_O,f_rate,r_rate,c_rate,wf_rate,wr_rate,matrix,res,react,dist);
						add_particle_info(particle,particle->get_type(),i,totals);
						//if (react[0]) fprintf(file4, "%i %f\n",particle->get_type(),dist);
						if (react[1] || react[2] || react[4]) delete particle;
						for (int k=0;k<8;k++) np_changes[k]+=res[k];
						f_reactions+=react[0];
						m++;
					}
					else {
						if (particle->get_interact() && (((double)rand()/(RAND_MAX))<=spin)) particle->spin_particle();
						add_particle_info(particle,particle->get_type(),i,totals);
					}
				}
			}
		}
		
		// Updates concentrations
		for (int k=0;k<8;k++) {
			nparticles[k]+=np_changes[k];
			if (k==4)
				conc_particles[k]=nparticles[k]*size_O/matrix_size;
			else
				conc_particles[k]=nparticles[k]/matrix_size;
		}
		
		gamma+=f_reactions/matrix_size;
		fprintf(file1,"%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",i,conc_particles[0],conc_particles[1],conc_particles[2],conc_particles[3],conc_particles[4],conc_particles[5],conc_particles[6],conc_particles[7],gamma);
		fprintf(file2,"%i\t %f\t %f\t %f\t %f\t %f\t %f\t",i,totals[1][0],totals[1][1],totals[1][2],totals[1][3],totals[1][4],totals[1][5]);
		fprintf(file2,"%f\t %f\t %f\t %f\t %f\t %f\t",totals[3][0],totals[3][1],totals[3][2],totals[3][3],totals[3][4],totals[3][5]);
		fprintf(file2,"%f\t %f\t %f\t %f\t %f\t %f\t\n",totals[4][0],totals[4][1],totals[4][2],totals[4][3],totals[4][4],totals[4][5]);
		total_nparticles=nparticles[0]+nparticles[1]+nparticles[2]+nparticles[3];
		if (mobil_O) total_nparticles+=nparticles[4];
		if (wfi>0) total_nparticles+=nparticles[5];
				
		//write_particles_pos(i,matrix,file3);
		i++;
	}
	
	fclose(file1); // Closes the results1 file
	fclose(file2); // Closes the results2 file
	//fclose(file3); // Closes the results3 file	
	//fclose(file4); // Closes the results4 file
	
	// Displays final results
	cout << nparticles[0] << " " << nparticles[1] << " " << nparticles[2] << " " << nparticles[3] << " " << nparticles[4] << " " << nparticles[5] << "\n";
	for (int x=1;x<5;x++) {
		for (int y=0;y<6;y++)
			cout << totals[x][y] << " ";
		cout << "\n";
	}
	cout << "\n\n";
	
    return 0;
}
