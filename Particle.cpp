/*
 *  Particle.cpp
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 2/5/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "Complex.h"
#include <stdlib.h>

// Constructors and Destructor
Particle::Particle(int new_id,int new_type,int new_pos_x,int new_pos_y,int new_it_creation){
	id=new_id;
	type=new_type;
	init_pos_x=new_pos_x;
	init_pos_y=new_pos_y;
	pos_x=new_pos_x;
	pos_y=new_pos_y;
	it_move=new_it_creation;
	it_creation=new_it_creation;
	part_move=0;
	sq_disp=0;
	pdist=0;
};

// Getters
int Particle::get_id() {
	return id;
};

int Particle::get_type() { 
	return type;
};

int Particle::get_pos_x() { 
	return pos_x;
};

int Particle::get_pos_y() {
	return pos_y;
};

int Particle::get_ori_x() { 
	return ori_x;
};

int Particle::get_ori_y() {
	return ori_y;
};

int Particle::get_ori_id() {
	return ori_id;
};

int Particle::get_init_pos_x() { 
	return init_pos_x;
};

int Particle::get_init_pos_y() {
	return init_pos_y;
};

int Particle::get_part_move() {
	return part_move;
};

double Particle::get_sq_disp() {
	return sq_disp;
};

int Particle::get_it_move() {
	return it_move;
};

int Particle::get_it_creation() {
	return it_creation;
};

double Particle::get_pdist() {
	return pdist;
};

bool Particle::get_interact() {
	return interact;
};

// Setters
void Particle::set_type(int new_type) {
	type=new_type;
}

void Particle::set_pos_x(int new_pos_x) {
	pos_x=new_pos_x;
}

void Particle::set_pos_y(int new_pos_y) {
	pos_y=new_pos_y;
}

void Particle::set_ori_x(int new_ori_x) {
	ori_x=new_ori_x;
}

void Particle::set_ori_y(int new_ori_y) {
	ori_y=new_ori_y;
}

void Particle::set_ori_id(int new_ori_id) {
	ori_id=new_ori_id;
}

void Particle::set_init_pos_x(int new_init_pos_x) {
	init_pos_x=new_init_pos_x;
}

void Particle::set_init_pos_y(int new_init_pos_y) {
	init_pos_y=new_init_pos_y;
}

void Particle::set_pos_x_y(int new_pos_x,int new_pos_y) {
	pos_x=new_pos_x;
	pos_y=new_pos_y;
}

void Particle::set_it_move(int new_it_move) {
	it_move=new_it_move;
}

void Particle::set_it_creation(int new_it_creation) {
	it_creation=new_it_creation;
}

void Particle::set_pdist(double new_pdist) {
	pdist=new_pdist;
}

void Particle::set_interact(bool new_interact) {
	interact=new_interact;
}

void Particle::add_part_move(int new_add_move) {
	part_move+=new_add_move;
}

void Particle::add_sq_disp(double new_add_sq_disp) {
	sq_disp+=new_add_sq_disp;
}

void Particle::create_ori_momentum() {
	ori_mom=pow(-1,rand()%2);
}

// Other Operations

void Particle::display_particle() {
	std::cout << "TYPE: " << type << " ";
	std::cout << "Pos_X: " << pos_x << "  ";
	std::cout << "Pos_Y: " << pos_y << "  ";
}

// this function initializes the particle's orientation
void Particle::init_spin_particle() {
	int new_ori_x[6]={1,1,0,-1,-1,0};
	int new_ori_y[6]={0,-1,-1,0,1,1};
	
	ori_id=rand()%6;
	ori_x=new_ori_x[ori_id];
	ori_y=new_ori_y[ori_id];
	create_ori_momentum();
}

// This function spins the particle
void Particle::spin_particle() {
	int new_ori_x[6]={1,1,0,-1,-1,0};
	int new_ori_y[6]={0,-1,-1,0,1,1};
	
	ori_id=ori_id+ori_mom;
	if (ori_id==-1) ori_id=5;
	if (ori_id==6) ori_id=0;
	ori_x=new_ori_x[ori_id];
	ori_y=new_ori_y[ori_id];
}

// This function checks the boundaries on the matrix
int Particle::check_boundaries(int final_axis_pos,int matrix_dim_size,bool& cross) {
	if (final_axis_pos<0) {
		final_axis_pos=matrix_dim_size-abs(final_axis_pos);
		cross=1;
	}
	else if (final_axis_pos>matrix_dim_size-1) {
		final_axis_pos=final_axis_pos-matrix_dim_size;
		cross=1;
	}
	return final_axis_pos;
}

// This function chooses a square destination for the particle (new function)
void Particle::choose_destination_square(int* move,int* res,bool& cross) {
	move[0]=0;
	move[1]=0;
	
	// Selects a destination where to move to
	int rand_value=rand()%4; // produces a random number between 0 and 3
	switch (rand_value) {
		case 0: {move[0]=-1;move[1]=0;break;}
		case 1: {move[0]=0;move[1]=1;break;}
		case 2: {move[0]=1;move[1]=0;break;}
		case 3: {move[0]=0;move[1]=-1;}
	}
	res[0]=this->check_boundaries(pos_x+move[0],matrix_height,cross);
	res[1]=this->check_boundaries(pos_y+move[1],matrix_width,cross);
}

// This function chooses a triangular destination for the particle (new function)
void Particle::choose_destination_triangular(double spin,int* move,int* res,bool& cross) {
	move[0]=0;
	move[1]=0;
	
	if (spin>0 && ((type==2) || type==5)) {
		move[0]=this->get_ori_x();
		move[1]=this->get_ori_y();
	}
	else {
		int rand_value=rand()%6; // produces a random number between 0 and 5
		switch (rand_value) {
			case 0: {move[0]=1;move[1]=0;break;}
			case 1: {move[0]=1;move[1]=-1;break;}
			case 2: {move[0]=0;move[1]=-1;break;}
			case 3: {move[0]=-1;move[1]=0;break;}
			case 4: {move[0]=-1;move[1]=1;break;}
			case 5: {move[0]=0;move[1]=1;}
		}
	}
	
	res[0]=this->check_boundaries(pos_x+move[0],matrix_height,cross);
	res[1]=this->check_boundaries(pos_y+move[1],matrix_width,cross);
}

// This function moves a particle from one spot of the grid to another
void Particle::move_particle(int* final_pos,double spin,bool cross,Particle* matrix[][matrix_width]) {
	if (cross) {
		pdist=pdist+sqrt(pow(final_pos[0]-init_pos_x,2)+pow(final_pos[1]-init_pos_y,2));
		pdist=pdist+sqrt(pow(final_pos[0]-init_pos_x,2)+pow(final_pos[1]-init_pos_y,2));
		init_pos_x=final_pos[0];
		init_pos_y=final_pos[1];
	}
	add_sq_disp(1);
	add_part_move(1);
	matrix[pos_x][pos_y]=NULL;
	pos_x=final_pos[0];
	pos_y=final_pos[1];
	matrix[final_pos[0]][final_pos[1]]=this;
	if (((double)rand()/(RAND_MAX))<=spin) spin_particle(); // Spins the particle to converse the momentum
}

// This function creates a complex
void Particle::create_complex(int new_type,int iteration,int* final_pos,Particle* dest_part,Particle* matrix[][matrix_width]){
	Complex* comp=NULL;
	Particle* part=NULL;
	
	if (((new_type==2) && (type==0)) || ((new_type==5) && (type==4))) {
		comp=new Complex(0,new_type,final_pos[0],final_pos[1],iteration,this,dest_part);
		part=this;
	}
	else {
		comp=new Complex(0,new_type,final_pos[0],final_pos[1],iteration,dest_part,this);
		part=dest_part;
	}
	
	comp->set_ori_x(part->get_ori_x());
	comp->set_ori_y(part->get_ori_y());
	comp->set_ori_id(part->get_ori_id());
	comp->create_ori_momentum();
	matrix[pos_x][pos_y]=NULL;
	matrix[final_pos[0]][final_pos[1]]=comp;
}

// This is the main function of class particle
void Particle::move(int iteration,double spin,double mobil_O,double f_rate,double r_rate,double c_rate,double wf_rate,double wr_rate,Particle* matrix[][matrix_width],int* res,int* react,double& dist) {
	bool cross=0;
	int new_type=0;
	int move[]={0,0};
	int final_pos[]={0,0};
	for (int k=0;k<8;k++) res[k]=0;
	for (int k=0;k<5;k++) react[k]=0;
	
	//choose_destination_square(move,final_pos,cross); // Determines the new destination using a triangular lattice
	choose_destination_triangular(spin,move,final_pos,cross); // Determines the new destination using a triangular lattice
	Particle* dest_part=matrix[final_pos[0]][final_pos[1]]; // Gets the particle at the destination site
	if (dest_part!=NULL) { // destination is occupied
		add_part_move(1);
		int dest_part_type=dest_part->get_type();
		
		// Checks if particles may react with each other and determines the type of the potential new particle
		if ((type<2) && (dest_part_type<2) && (type!=dest_part_type))
			new_type=2;
		else if ((type==1 && dest_part_type==4 && dest_part->get_interact()) || (type==4 && interact && dest_part_type==1)) // Original: type < 4 || dest_part_type < 4
			new_type=5;
		
		// Creates a new particle if all the conditions are met
		if (new_type!=0) {
			bool mismatch_ori=0;
			if (spin>0) mismatch_ori=ori_x+dest_part->get_ori_x()+ori_y+dest_part->get_ori_y();
			if (mismatch_ori==0) {
				double rn=((double)rand()/(double)RAND_MAX); // Generate a double random number between 0 and 1
				if ((new_type==2) && (rn<f_rate)){ // forward reaction takes place
					react[0]=1;
					add_sq_disp(1);
					res[0]=-1;res[1]=-1;res[2]=1;
					create_complex(new_type,iteration,final_pos,dest_part,matrix);
					dist=pdist+sqrt(pow(final_pos[0]-init_pos_x,2)+pow(final_pos[1]-init_pos_y,2));
				}
				else if ((new_type==5) && (rn<wf_rate)){ // forward reaction takes place
					add_sq_disp(1);
					react[3]=1;
					res[4]=-1;
					res[5]=1;
					if (type<4){
						res[type]=-1;
						if (type==0) res[6]=1;
						if (type==3) res[7]=1;
					}
					else {
						res[dest_part_type]=-1;
						if (dest_part_type==0) res[6]=1;
						if (dest_part_type==3) res[7]=1;
					}
					create_complex(new_type,iteration,final_pos,dest_part,matrix);
				}
			}
		}
		
		// Shifts the momentum and probably spins the particle if it collided with anything that did not react with
		if (!(react[0] || react[3])) {
			create_ori_momentum();
			if (((double)rand()/(RAND_MAX))<=spin) spin_particle();
			if (dest_part_type!=4 || dest_part->get_interact()) {
				dest_part->create_ori_momentum();
				if (((double)rand()/(RAND_MAX))<=spin) dest_part->spin_particle();
			}
		}
	}
	else if (type==2) { // Particle is a Complex (type 2)
		double rn=((double)rand()/(double)RAND_MAX); // Generate a double random number between 0 and 1
		if (rn<r_rate) { // reverse reaction takes place
			int rel_type;
			react[1]=1;
			res[0]=1;res[1]=1;res[2]=-1;
			reverse_reaction(final_pos,iteration,rel_type,matrix);
		}
		else if (rn<(r_rate+c_rate)) { // catalytic reaction takes place
			react[2]=1;
			res[0]=1;res[2]=-1;res[3]=1;
			catalytic_reaction(final_pos,iteration,matrix);
		}
		else move_particle(final_pos,spin,cross,matrix);
	}
	else if (type==5) { // Particle is a Complex (type 5)
		double rn=((double)rand()/(double)RAND_MAX); // Generate a double random number between 0 and 1
		if (rn<wr_rate) { // reverse reaction takes place
			int rel_type;
			react[4]=1;
			reverse_reaction(final_pos,iteration,rel_type,matrix);
			res[4]=1;
			res[5]=-1;
			res[rel_type]=1;
			if (rel_type==0) res[6]=-1;
			if (rel_type==3) res[7]=-1;
		}
		else if (mobil_O) move_particle(final_pos,spin,cross,matrix);
	}
	else move_particle(final_pos,spin,cross,matrix);
}
