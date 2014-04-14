/*
 * sphdem.cpp
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 7 Feb 2014
 *      Author: robinsonm
 */

#include "crowding.h"
#include "Visualisation.h"

#include <random>

int main(int argc, char **argv) {
	auto A = SpeciesType::New();
	auto params = ptr<Params>(new Params());


	const int timesteps = 20000;
	const int nout = 1000;
	const int timesteps_per_out = timesteps/nout;
	const int n = 100;

	/*

	/*
	 * parameters
	 */

	const double kDa = 1.660538921e-30;
	const double mass = 40.0*kDa; //average mol mass of protein in E.Coli
	const double viscosity = 8.9e-4; //viscosity of water
	params->diameter = 5e-9; //average diameter of protein in E.Coli
	params->time = 0;
	params->T = 300.0; //room temp
	params->D = k_b*params->T/(3.0*PI*viscosity*params->diameter);
	params->k_s = 1.0;
	const double aim_step_length = params->diameter/100.0;
	params->dt =  pow(aim_step_length,2)/(2.0*params->D);
	const double gamma = 3.0*PI*viscosity*params->diameter/mass;


	const double L = params->diameter*10;



	std::cout <<"Running simulation with parameters:"<<std::endl;
	std::cout <<"\tD = "<<params->D<<" average step length = "<<sqrt(2.0*params->D*params->dt)<<" asl = "<<sqrt(2.0*params->D*params->dt)/L<<" L"<<std::endl;
	std::cout <<"\tdiameter = "<<params->diameter<<std::endl;
	std::cout <<"\tdt = "<<params->dt<<" m/gamma = "<<mass/gamma<<" mean free time = "<<params->diameter/(2.0*sqrt(3.0*k_b*params->T/mass))<<std::endl;
	std::cout <<"\tT = "<<params->T<<std::endl;
	std::cout <<"\tk_s = "<<params->k_s<<std::endl;


	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,true);

	/*
	 * create particles
	 */
	std::mt19937 generator;
	std::uniform_real_distribution<double> distribution(0,L);
	auto dice = std::bind ( distribution, generator );
	A->create_particles(n,[L,A,&dice](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		v << 0,0,0;
		Vect3d canditate_position;
		bool regenerate;
		do {
			regenerate = false;
			canditate_position << dice(),dice(),dice();
			for (auto tpl: A->get_neighbours(canditate_position)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				if (alivej) {
					regenerate = true;
					break;
				}
			}
		} while (regenerate);
		return canditate_position;
	});


	/*
	 * setup output stuff
	 */
	auto A_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	A->copy_to_vtk_grid(A_grid);
	Visualisation::vtkWriteGrid("vis/at_start_A",0,A_grid);

	A->init_neighbour_search(min,max,params->diameter,periodic);

	/*
	 * Simulate!!!!
	 */
	std::cout << "starting...."<<std::endl;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			langevin_timestep(A,params);
		}
		std::cout <<"iteration "<<i<<std::endl;
		
		A->copy_to_vtk_grid(A_grid);
		Visualisation::vtkWriteGrid("vis/A",i,A_grid);
	}
	
	
}
