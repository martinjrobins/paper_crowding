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


	if (argc < 8) {
		std::cout << "Usage: self_crowding output_dir time nout k_s aim_step_length/diameter vol_ratio method"<<std::endl;
		return -1;
	}


	const double sim_time = atof(argv[2]);
	const int nout = atoi(argv[3]);

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
	//params->k_s = 1.0e-2;
	params->k_s = atof(argv[4]);
	const double aim_step_length = params->diameter*atof(argv[5]);
	params->dt =  pow(aim_step_length,2)/(2.0*params->D);
	const int timesteps = sim_time/params->dt;
	const int timesteps_per_out = timesteps/nout;
	const double gamma = 3.0*PI*viscosity*params->diameter/mass;


	const double L = params->diameter*20;
	const double rdf_min = params->diameter*0.1;
	const double rdf_max = params->diameter*3;
	const int rdf_n = 100;

	const double vol_ratio = atof(argv[6]);
	const double mol_vol = (1.0/6.0)*PI*pow(params->diameter,3);
	const int n = pow(L,3)*vol_ratio/mol_vol;



	std::cout <<"Running simulation with parameters:"<<std::endl;
	std::cout <<"\tvolume ratio = "<<vol_ratio<<" n = "<<n<<std::endl;
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
	A->create_particles(n,[params,L,A,&dice](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		Vect3d canditate_position;
		bool regenerate;
		do {
			regenerate = false;
			canditate_position << dice(),dice(),dice();
			for (auto tpl: A->get_neighbours(canditate_position)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > params->diameter*params->diameter) continue;
				if (alivej) {
					regenerate = true;
					break;
				}
			}
		} while (regenerate);
		r0 = canditate_position;
		rt = canditate_position;
		v << 0,0,0;
		U = 0;
		exits = 0;
		return canditate_position;
	});


	/*
	 * setup output stuff
	 */
	auto A_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	A->copy_to_vtk_grid(A_grid);
	std::string output_dir(argv[1]);
	Visualisation::vtkWriteGrid((output_dir+"/at_start_A").c_str(),0,A_grid);

	A->init_neighbour_search(min,max,params->diameter,periodic);

	std::vector<double> rdf_r;
	rdf_r.resize(rdf_n);
	for (int i = 0; i < rdf_n; ++i) {
		rdf_r[i] = i*(rdf_max-rdf_min)/rdf_n + rdf_min;
	}

	std::ofstream f;
	f.open((output_dir+"/msd.csv").c_str());
	f << "#timestep,time,msd,flux"<<std::endl;

	enum method{LANGEVIN,MONTE_CARLO,HARD_SPHERE};
	enum method the_method = LANGEVIN;
	if (strcmp(argv[7],"langevin")==0) the_method=LANGEVIN;
	if (strcmp(argv[7],"monte_carlo")==0) the_method=MONTE_CARLO;
	if (strcmp(argv[7],"hard_sphere")==0) the_method=HARD_SPHERE;
	std::cout << the_method<<std::endl;

	/*
	 * Simulate!!!!
	 */
	std::cout << "starting...."<<std::endl;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			switch (the_method) {
			case LANGEVIN:
				langevin_timestep(A,params);
				break;
			case MONTE_CARLO:
				monte_carlo_timestep(A,params,generator);
				break;
			case HARD_SPHERE:
				hard_sphere_timestep(A,params);
				break;
			default:
				break;
			}
		}
		std::cout <<"iteration "<<i<<std::endl;
		A->copy_to_vtk_grid(A_grid);
		Visualisation::vtkWriteGrid((output_dir+"/A").c_str(),i,A_grid);
		if (i==10) {
			std::for_each(A->begin(),A->end(),[](SpeciesType::Value& i) {
				REGISTER_SPECIES_PARTICLE(i);
				r0 = r;
				rt = r;
				exits = 0;
			});
		}
		double msv = std::accumulate(A->begin(),A->end(),0.0,[](double i,SpeciesType::Value& j) {
			const GET_TUPLE(Vect3d,rtj,SPECIES_TOTAL_R,j);
			const GET_TUPLE(Vect3d,r0j,SPECIES_SAVED_R,j);
			return i + (rtj-r0j).squaredNorm();
		})/A->size();

		const double flux = std::accumulate(A->begin(),A->end(),0.0,[](unsigned int i,SpeciesType::Value& j) {
			const GET_TUPLE(unsigned int,exits,SPECIES_NUM_EXITS,j);
			return i + exits;
		})/(pow(L,2)*6.0);

		f << i<<","<<(i+1)*timesteps_per_out*params->dt<<","<<msv<<','<<flux<<std::endl;

		auto rdf = radial_distribution_function(A,rdf_min,rdf_max,rdf_n);
		char buffer[100];
		sprintf(buffer,"%s/rdf%05d.csv",argv[1],i);
		Visualisation::write_column_vectors(buffer,"#r,rdf",{rdf_r,*rdf});

	}
	f.close();

	
}
