/*
 *  Particle.h
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 2/5/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <math.h>

#define matrix_width 100 // sets the width of the matrix
#define matrix_height 100 // sets the height of the matrix

class Particle {
	
public:
	
	// Constructors and Destructors
	Particle(int,int,int,int,int);
	~Particle(){};
	
	// Getters
	int get_id();
	int get_type();
	int get_pos_x();
	int get_pos_y();
	int get_ori_x();
	int get_ori_y();
	int get_ori_id();
	int get_init_pos_x();
	int get_init_pos_y();
	int get_part_move();
	double get_sq_disp();
	bool get_interact();
	int get_it_move();
	int get_it_creation();
	double get_pdist();
	
	// Setters
	void set_type(int);
	void set_pos_x(int);
	void set_pos_y(int);
	void set_ori_x(int);
	void set_ori_y(int);
	void set_ori_id(int);
	void set_init_pos_x(int);
	void set_init_pos_y(int);
	void set_pos_x_y(int,int);
	void set_it_move(int);
	void set_it_creation(int);
	void set_pdist(double);
	void set_interact(bool);
	void add_part_move(int);
	void add_sq_disp(double);
	void create_ori_momentum();
	
	// Other Operations
	void spin_particle();
	void display_particle();
	void init_spin_particle();
	int check_boundaries(int,int,bool&); // Used in Obstacle
	void move_particle(int*,double,bool,Particle* [][matrix_width]);
	void choose_destination_square(int*,int*,bool&); // Used in Obstacle
	void choose_destination_triangular(double,int*,int*,bool&); // Used in Obstacle
	void create_complex(int,int,int*,Particle*,Particle*[][matrix_width]);
	virtual void reverse_reaction(int*,int,int&,Particle* [][matrix_width]){}; // Defined in Complex
	virtual void catalytic_reaction(int*,int,Particle* [][matrix_width]){}; // Defined in Complex
	void move(int,double,double,double,double,double,double,double,Particle* [][matrix_width],int*,int*,double&);

protected:
	
	int id; // Stores the particle's id
	int type; // Stores the particle's type
	int pos_x; // Stores position of the particle in the x axis
	int pos_y; // Stores position of the particle in the y axis
	int ori_x; // stores the particle's orientation x position
	int ori_y; // stores the particle's orientation y position
	int ori_id; // Stores the particle's orientation id
	int ori_mom; // Stores the particle's orientational momentum
	int part_move; // Stores the number of times the particle has moved
	double sq_disp; // Stores the square displacement for the particles

private:
	
	int init_pos_x; // Stores initial position of the particle in the x axis
	int init_pos_y; // Stores initial position of the particle in the x axis
	int it_move; // Stores the iteration at which the particle has moved
	int it_creation; // Stores the iterations at which the particles were created
	double pdist; // Stores the distance travelled by the particle
	bool interact; // States the obstacle particle is interacting
};

