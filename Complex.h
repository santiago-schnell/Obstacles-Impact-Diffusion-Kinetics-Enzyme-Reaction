/*
 *  Complex.h
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 3/5/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "Particle.h"

class Complex : public Particle {
	
public:
	// Constructors and Destructors
	Complex(int,int,int,int,int,Particle*,Particle*);

	void initialize_new_particle(Particle*,int,int,int);
	void reverse_reaction(int*,int,int&,Particle* [][matrix_width]);
	void catalytic_reaction(int*,int,Particle* [][matrix_width]);
	
private:
	Particle* assoc_particles[2]; // Stores the complex particles
};