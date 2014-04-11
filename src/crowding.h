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

#ifndef SPHDEM_H_
#define SPHDEM_H_

#include "Aboria.h"
using namespace Aboria;

#include <tuple>

const double k_b = 1.3806488e-23;

enum {SPECIES_VELOCITY};
typedef std::tuple<Vect3d> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPECIES_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				const bool alive = particle.is_alive(); \
				GET_TUPLE(Vect3d,v,SPECIES_VELOCITY,particle);
#define REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const DemType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const bool alivej = j.is_alive(); \
				const GET_TUPLE(Vect3d,vj,SPECIES_VELOCITY,j);

struct Params {
	double diameter,D,dt,T,k_s;
	double time;
};


template<typename GeometryType>
void timestep(ptr<SpeciesType> A,
		ptr<Params> params) {

	const double dt = params->dt;
	const double diameter = params->diameter;
	const double D = params->D;
	const double T = params->T;
	const double k_s = params->k_s;

	A->update_positions(A->begin(),A->end(),[dt,diameter,D,T,k_b,k_s](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		v << 0,0,0;
		for (auto tpl: i.get_neighbours(A)) {
			REGISTER_NEIGHBOUR_SPECES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (r2 == 0) continue;

			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				v += (k_s*overlap)*normal;
			}
		}
		v *= dt*D/(k_b*T);
		v += sqrt(2.0*D*dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());

		return r + v;
	});

}


#endif /* SPHDEM_H_ */
