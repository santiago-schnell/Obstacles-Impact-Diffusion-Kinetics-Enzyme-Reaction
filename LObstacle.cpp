/*
 *  LObstacle.cpp
 *  Diffusion2
 *
 *  Created by Marcio Duarte Albasini Mourao on 3/6/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "LObstacle.h"

LObstacle::LObstacle(int a,int b,int c,int d,int e,int f,int g):Particle(a,b,c,d,e) {
	init_obst_pos_x=f;
	init_obst_pos_y=g;
}

// Getters
int LObstacle::get_init_obst_pos_x() { 
	return init_obst_pos_x;
};

int LObstacle::get_init_obst_pos_y() { 
	return init_obst_pos_y;
};

// Setters
void LObstacle::set_init_obst_pos_x(int new_init_obst_pos_x) {
	init_obst_pos_x=new_init_obst_pos_x;
}

void LObstacle::set_init_obst_pos_y(int new_init_obst_pos_y) {
	init_obst_pos_y=new_init_obst_pos_y;
}

// This function checks the large obstacle for creation errors
bool check_obstacle(LObstacle* this_part,int size,int init_x,int init_y,Particle* matrix[][matrix_width]) {
	Particle* part=NULL;
	int id=0,j=0,k=0;
	int posx,posy;
	bool errors=0;
	bool cross=0;
	
	part=matrix[init_x][init_y];
	if ((part!=NULL) && (part->get_type()==4)) {
		id=part->get_id();
		j=0;
		while ((errors==0) && (j<sqrt(size))) {
			k=0;
			while ((errors==0) && (k<sqrt(size))) {
				posx=this_part->check_boundaries(init_x+j,matrix_height,cross);
				posy=this_part->check_boundaries(init_y+k,matrix_width,cross);
				if ((matrix[posx][posy]==NULL) || (matrix[posx][posy]->get_id()!=id))
					errors=1;
				k++;
			}
			j++;
		}
	}
	return errors;
}

// This function moves an obstacle of n size
void move_obstacle(LObstacle* this_part,int init_obst_pos_x,int init_obst_pos_y,int* move,int size,bool& cross,Particle* matrix[][matrix_width]) {
	LObstacle* particles[2500]={NULL};
	int count=0;
	int posx,posy;
	int old_posx,old_posy;
	int new_init_obst_pos_x,new_init_obst_pos_y;
	int hlp1,hlp2;
	
	hlp1=move[0];
	hlp2=move[1];
	
	for (int j=0;j<sqrt(size);j++)
		for (int k=0;k<sqrt(size);k++) {
			old_posx=this_part->check_boundaries(init_obst_pos_x+j,matrix_height,cross);
			old_posy=this_part->check_boundaries(init_obst_pos_y+k,matrix_width,cross);
			particles[count]=(LObstacle*)matrix[old_posx][old_posy];
			posx=this_part->check_boundaries(init_obst_pos_x+j+move[0],matrix_height,cross);
			posy=this_part->check_boundaries(init_obst_pos_y+k+move[1],matrix_width,cross);
			if ((j==0) && (k==0)) {
				new_init_obst_pos_x=posx;
				new_init_obst_pos_y=posy;
			}
			particles[count]->set_init_obst_pos_x(new_init_obst_pos_x);
			particles[count]->set_init_obst_pos_y(new_init_obst_pos_y);
			particles[count]->set_pos_x(posx);
			particles[count]->set_pos_y(posy);
			particles[count]->add_part_move(1);
			particles[count]->add_sq_disp(1);
			matrix[old_posx][old_posy]=NULL;
			count++;
		}
	
	for (int j=0;j<count;j++) {
		posx=particles[j]->get_pos_x();
		posy=particles[j]->get_pos_y();
		matrix[posx][posy]=particles[j];
	}
}			

// this function communicates the intention to move
void add_part_move_obst_particles(LObstacle* this_part,int size,int init_obst_pos_x,int init_obst_pos_y,Particle* matrix[][matrix_width]) {
	bool cross=0;
	int posx,posy;
	
	for (int j=0;j<sqrt(size);j++)
		for (int k=0;k<sqrt(size);k++) {
			posx=this_part->check_boundaries(init_obst_pos_x+j,matrix_height,cross);
			posy=this_part->check_boundaries(init_obst_pos_y+k,matrix_width,cross);
			matrix[posx][posy]->add_part_move(1);
		}
}

// this function checks if the large obstacle can move
bool check_move_availability(Particle* this_part,int id,int init_obst_pos_x,int init_obst_pos_y,int* move,int size,bool& cross,Particle* matrix[][matrix_width]) {
	bool avail=1;
	int j=0,k=0;
	int posx,posy;
	Particle* particle;
	
	while (avail && (j<sqrt(size))) {
		k=0;
		while (avail && (k<sqrt(size))) {
			posx=this_part->check_boundaries(init_obst_pos_x+j+move[0],matrix_height,cross);
			posy=this_part->check_boundaries(init_obst_pos_y+k+move[1],matrix_width,cross);
			particle=matrix[posx][posy];
			if ((particle!=NULL) && (particle->get_id()!=id))
				avail=0;
			k++;
		}
		j++;
	}
	return avail;
}

// This is the main function of class LObstacle
void LObstacle::move(int size,Particle* matrix[][matrix_width]) {
	bool res=1;
	bool cross=0;
	int move[]={0,0};
	int final_pos[]={0,0};
		
	//choose_destination_square(move,final_pos,cross);
	choose_destination_triangular(0,move,final_pos,cross);
	res=check_move_availability(this,id,init_obst_pos_x,init_obst_pos_y,move,size,cross,matrix);
	if (res)
		move_obstacle(this,init_obst_pos_x,init_obst_pos_y,move,size,cross,matrix);
	else
		add_part_move_obst_particles(this,size,init_obst_pos_x,init_obst_pos_y,matrix);
}
