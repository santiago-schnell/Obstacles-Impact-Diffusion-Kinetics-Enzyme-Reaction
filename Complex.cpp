/*
 *  Complex.cpp
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 3/5/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "Complex.h"

Complex::Complex(int a,int b,int c,int d,int e,Particle* part1,Particle* part2):Particle(a,b,c,d,e){
	assoc_particles[0]=part1;
	assoc_particles[1]=part2;
};

void Complex::initialize_new_particle(Particle* particle,int new_x,int new_y,int iteration) {
	int new_type=particle->get_type();
	particle->set_pos_x(new_x);
	particle->set_pos_y(new_y);
	particle->create_ori_momentum();

	if ((type==5) || (new_type!=3)) {
		particle->set_it_move(iteration);
		particle->add_part_move(part_move);
		particle->add_sq_disp(sq_disp);
	}
	
	if (((type==2) && (new_type==0)) || (new_type==4)){
		particle->set_ori_x(ori_x);
		particle->set_ori_y(ori_y);
		particle->set_ori_id(ori_id);
	}
	else {
		particle->set_ori_x(0);
		particle->set_ori_y(0);
		if (ori_x!=0) particle->set_ori_x((-1)*ori_x);
		if (ori_y!=0) particle->set_ori_y((-1)*ori_y);
		if (ori_id<3)
			particle->set_ori_id(ori_id+3);
		else
			particle->set_ori_id(ori_id-3);
	}
}

void Complex::reverse_reaction(int* final_pos,int iteration,int& rel_type,Particle* matrix[][matrix_width]){
	Particle* part1=assoc_particles[0]; // This particle will stay in place
	Particle* part2=assoc_particles[1]; // This particle will move
	
	initialize_new_particle(part1,pos_x,pos_y,iteration);
	initialize_new_particle(part2,final_pos[0],final_pos[1],iteration);
	rel_type=assoc_particles[1]->get_type();
	part2->add_part_move(1);
	part2->add_sq_disp(1);
	matrix[pos_x][pos_y]=part1;
	matrix[final_pos[0]][final_pos[1]]=part2;
};

void Complex::catalytic_reaction(int* final_pos,int iteration,Particle* matrix[][matrix_width]){
	Particle* enzyme=assoc_particles[0];
	Particle* product = new Particle(0,3,final_pos[0],final_pos[1],iteration); // create a particle at the specified position
	
	initialize_new_particle(enzyme,pos_x,pos_y,iteration);
	initialize_new_particle(product,final_pos[0],final_pos[1],iteration);
	product->add_part_move(1);
	product->add_sq_disp(1);
	matrix[pos_x][pos_y]=enzyme;
	matrix[final_pos[0]][final_pos[1]]=product;
};
