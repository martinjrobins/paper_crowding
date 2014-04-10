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

#include <vtkFloatArray.h>


int main(int argc, char **argv) {
	auto A = SpeciesType::New();


	const int timesteps = 20000;
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 0.004;
	const int nx = 10;



	/*

	/*
	 * parameters
	 */



	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,true);

	/*
	 * create particles
	 */
	dem->create_particles(min,max,dem_n,[](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		v << 0,0,0;
		v0 << 0,0,0;
		f << 0,0,0;
		f0 << 0,0,0;
	});

	sph->create_particles_grid(min_domain,max,sph_n+Vect3i(0,0,3),[psep,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		h = params->sph_maxh;
		omega = 1.0;
		kappa = 0.0;
		v << 0,0,0;
		v0 << 0,0,0;
		dddt = 0;
		f0 << 0,0,0;
		e = 1;
		rho = params->sph_dens;
		f << 0,0,0;
		fdrag << 0,0,0;
		f0 << 0,0,0;
		if (r[2]<0) {
			fixed = true;
		} else {
			fixed = false;
		}
	});

	/*
	 * setup output stuff
	 */
	auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph->copy_to_vtk_grid(sph_grid);
	dem->copy_to_vtk_grid(dem_grid);
	Visualisation::vtkWriteGrid("vis/at_start_sph",0,sph_grid);
	Visualisation::vtkWriteGrid("vis/at_start_dem",0,dem_grid);


	std::cout << "init porosity and rho...."<<std::endl;
	sph->init_neighbour_search(min_domain,max_domain,2*params->sph_maxh,periodic);
	dem->init_neighbour_search(min_domain,max_domain,2*params->sph_maxh,periodic);

	/*
	 * init porosity and rho
	 */
	//std::cout << "calculate omega"<<std::endl;
	std::for_each(sph->begin(),sph->end(),[sph,dem,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		e = 1;
		//bool found = false;
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			const double r = sqrt(r2);
			const double q = r/h;
			const double Wab = W(q,h);
			e -= params->dem_vol*Wab;
			//found = true;
		}
		rho = params->sph_dens/e;
	});


	/*
	 * Simulate!!!!
	 */
	std::cout << "starting...."<<std::endl;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			sphdem(sph,dem,params,sph_geometry,dem_geometry);
		}
		std::cout <<"iteration "<<i<<std::endl;
		
		sph->copy_to_vtk_grid(sph_grid);
		dem->copy_to_vtk_grid(dem_grid);
		Visualisation::vtkWriteGrid("vis/sph",i,sph_grid);
		Visualisation::vtkWriteGrid("vis/dem",i,dem_grid);
	}
	
	
}
