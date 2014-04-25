/*
 * sphdem.h
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

#ifndef CROWDING_H_
#define CROWDING_H_

#include "Aboria.h"
using namespace Aboria;

#include <tuple>

const double k_b = 1.3806488e-23;

enum {SPECIES_VELOCITY,SPECIES_POTENTIAL,SPECIES_TOTAL_R,SPECIES_SAVED_R,SPECIES_SAVED_R1,SPECIES_NUM_EXITS};
typedef std::tuple<Vect3d,double,Vect3d,Vect3d,Vect3d,unsigned int> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPECIES_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				const bool alive = particle.is_alive(); \
				GET_TUPLE(double,U,SPECIES_POTENTIAL,particle); \
				GET_TUPLE(unsigned int,exits,SPECIES_NUM_EXITS,particle); \
				GET_TUPLE(Vect3d,r0,SPECIES_SAVED_R,particle); \
				GET_TUPLE(Vect3d,rt,SPECIES_TOTAL_R,particle); \
				GET_TUPLE(Vect3d,r1,SPECIES_SAVED_R1,particle); \
				GET_TUPLE(Vect3d,v,SPECIES_VELOCITY,particle);
#define REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const SpeciesType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const bool alivej = j.is_alive(); \
				const GET_TUPLE(double,Uj,SPECIES_POTENTIAL,j); \
				const GET_TUPLE(Vect3d,vj,SPECIES_VELOCITY,j);

struct Params {
	double diameter,D,dt,T,k_s;
	double time;
};


void langevin_timestep(ptr<SpeciesType> A,
		ptr<Params> params) {

	const double dt = params->dt;
	const double diameter = params->diameter;
	const double D = params->D;
	const double T = params->T;
	const double k_s = params->k_s;

	A->update_positions(A->begin(),A->end(),[A,dt,diameter,D,T,k_s](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		v << 0,0,0;
		U = 0;
		for (auto tpl: i.get_neighbours(A)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (r2 == 0) continue;

			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				v += (k_s*overlap)*normal;
				U += 0.5*k_s*pow(overlap,2);
			}
		}
		v *= dt*D/(k_b*T);
		v += sqrt(2.0*D*dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());
		rt += v;
		const Vect3d new_position = r+v;
		if ((new_position.array() < A->get_low().array()).any() ||
				(new_position.array() >= A->get_high().array()).any()) {
			exits++;
		}
		return new_position;
	});

}

double calculate_total_potential(ptr<SpeciesType> A,
		ptr<Params> params) {
	double sumU;
	const double dt = params->dt;
	const double diameter = params->diameter;
	const double D = params->D;
	const double T = params->T;
	const double k_s = params->k_s;

	/*
	 * calculate potential
	 */
	std::for_each(A->begin(),A->end(),[A,dt,diameter,D,T,k_s](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		U = 0;
		for (auto tpl: i.get_neighbours(A)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (r2 == 0) continue;

			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				U += 0.5*k_s*pow(overlap,2);
			}
		}
	});

	sumU = std::accumulate(A->begin(),A->end(), 0.0,[](double i,SpeciesType::Value& j) {
		const GET_TUPLE(double,Uj,SPECIES_POTENTIAL,j);
		return i+Uj;
	});
	return sumU;
}

template<typename RT>
void monte_carlo_timestep(ptr<SpeciesType> A,
		ptr<Params> params, RT& generator) {

	const double dt = params->dt;
	const double diameter = params->diameter;
	const double D = params->D;
	const double T = params->T;
	const double k_s = params->k_s;

	std::uniform_real_distribution<double> uniformd(0,1);
	std::normal_distribution<double> normald(0,1);
	const Vect3d low = A->get_low();
	const Vect3d high = A->get_high();
	int index;
	for (int ii = 0; ii < A->size(); ++ii) {
	while(1) {
		/*
		 * generate new state x'
		 */
		const int index = uniformd(generator)*A->size();
		REGISTER_SPECIES_PARTICLE(((*A)[index]));
		Vect3d canditate_pos = r+sqrt(2.0*D*dt)*Vect3d(normald(generator),normald(generator),normald(generator));
		for (int d = 0; d < 3; ++d) {
			while (canditate_pos[d]<low[d]) {
				canditate_pos[d] += (high[d]-low[d]);
			}
			while (canditate_pos[d]>=high[d]) {
				canditate_pos[d] -= (high[d]-low[d]);
			}
		}
		double Udiff = 0;
		for (auto tpl: A->get_neighbours(r)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (j.get_id()==((*A)[index]).get_id()) continue;
			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				Udiff -= k_s*pow(overlap,2);
			}
		}
		for (auto tpl: A->get_neighbours(canditate_pos)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (j.get_id()==((*A)[index]).get_id()) continue;
			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				Udiff += k_s*pow(overlap,2);
			}
		}

		const double acceptance_ratio = exp(-Udiff/(k_b*T));
		//std::cout <<"dU = "<<newU-oldU<<" acceptance_ratio = "<<acceptance_ratio<<std::endl;
		if (uniformd(generator)<acceptance_ratio) {
			//std::cout <<"accepted"<<std::endl;

			A->update_positions(A->begin()+index,A->begin()+index+1,[canditate_pos](SpeciesType::Value& i) {
				return canditate_pos;
			});
			break;
		}
	}
	}
}


#endif /* CROWDING_H_ */
