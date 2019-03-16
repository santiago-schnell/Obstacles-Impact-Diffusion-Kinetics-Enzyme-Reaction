/*
 *  LObstacle.h
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 3/6/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "Particle.h"

class LObstacle : public Particle {
	
public:
	
	// Constructors and Destructors
	LObstacle(int,int,int,int,int,int,int);

	// Getters
	int get_init_obst_pos_x();
	int get_init_obst_pos_y();
	
	// Setters
	void set_init_obst_pos_x(int);
	void set_init_obst_pos_y(int);
	
	// Others
	void move(int size,Particle* matrix[][matrix_width]);
	
private:
	
	// Properties
	int init_obst_pos_x; // Stores initial position of the large obstacle particle in the x axis
	int init_obst_pos_y; // Stores initial position of the large obstacle particle in the y axis	
};