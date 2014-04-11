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
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 1.0;
	const int n = 10;

	/*

	/*
	 * parameters
	 */
	params->D = 1.0;
	params->diameter = 0.01;
	params->dt = 0.0001;
	params->time = 0;

	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,true);

	/*
	 * create particles
	 */
	std::mt19937 generator;
	std::uniform_real_distribution distribution(0,L);
	auto dice = std::bind ( distribution, generator );
	A->create_particles(n,[A,&dice](SpeciesType::Value& i) {
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
			timestep(A,params);
		}
		std::cout <<"iteration "<<i<<std::endl;
		
		A->copy_to_vtk_grid(A_grid);
		Visualisation::vtkWriteGrid("vis/A",i,A_grid);
	}
	
	
}
