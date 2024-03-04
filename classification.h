/*  Copyright (C) 2021 by Diego Conti, diego.conti@unimib.it      
                                                                     
    This file is part of Gleipnir
	                                                                     
    Gleipnir is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gleipnir is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gleipnir.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H

#include <wedge/wedge.h>
#include <filesystem>
using namespace Wedge;
using std::string;

/** Represents an immutable list of mathematical objects corresponding to a classification */
template <typename Classified> class Classification {
	vector<unique_ptr<Classified>> elements;
	vector<string> names;
protected:
/** Add an entry to the classification

 @param element The entry to be added. Ownership is transferred to this container
*/
	void add (unique_ptr<Classified>&& element) {elements.emplace_back(move(element)); names.push_back(ToString(elements.size()));}
/** Add an entry to the classification

 @param element The entry to be added. Ownership is transferred to this container
 @param name A name for this entry
*/
	void add (unique_ptr<Classified>&& element, string name) {elements.push_back(move(element)); names.push_back(name);}
	Classification() {}
public:
	template<typename T> Classification(std::initializer_list<T> args_for_elems) {
		for (auto& args: args_for_elems)
			elements.emplace_back(args);
	}

	const Classified& entry (OneBased index) {return *elements[index-1];}
	const string& name (OneBased index) {return names[index-1];}
	auto begin() const {return elements.begin();}
	auto end() const {return elements.end();}
};


template<typename Parameter, typename... Name>
lst make_symbols(Name&&... name) {
	return lst{Parameter{std::forward<Name>(name)}...};
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(LieAlgebraParameter,realsymbol)

lst symbols() {
	lst symbols=make_symbols<LieAlgebraParameter>(
		N.a,N.b,N.c,N.d,N.e,N.f,N.g,N.h,N.i,N.j,N.k,N.l,N.m,N.n,N.o,N.p,N.q,N.r,N.s,N.t,N.u,N.v,N.w,N.x,N.y,N.z,
		N.A,N.B,N.C,N.D,N.E,N.F,N.G,N.H,N.I,N.J,N.K,N.L,N.M,N.N,N.O,N.P,N.Q,N.R,N.S,N.T,N.U,N.V,N.W,N.X,N.Y,N.Z,
		N.alpha,N.beta,N.gamma,N.Gamma,N. delta,N.Delta,N.epsilon,N.zeta,N.eta,N.theta,N.Theta,N.kappa,N.lambda,N.Lambda,
		N.mu,N.nu,N.xi,N.Xi,N.rho,N.pi,N.Pi,N.sigma,N.Sigma,N.tau,N.upsilon,N.Upsilon,N.phi,N.Phi,N.chi,N.psi,N.Psi,N.omega,N.Omega
	);
	for (int i=0;i<20;++i)
		symbols.append(LieAlgebraParameter(N.a(i)));
	return symbols;
}	

class GroupClassification : public Classification<LieGroup> {
protected:
	using Classification<LieGroup>::add;
	void add(const char* structure_constants) {
		Classification<LieGroup>::add(make_unique<AbstractLieGroup<false>>	(structure_constants));
	}
	template <typename... T>
	void add(const char* structure_constants, const T&... names) {
		Classification<LieGroup>::add(make_unique<AbstractLieGroup<true>>	(structure_constants, names...));
	}
public:
	static GroupClassification from_disk(const string& filename) {
		std::filesystem::path p{filename};		
		if (!std::filesystem::exists(p))
			p=std::filesystem::current_path().parent_path()/filename;
		if (!std::filesystem::exists(p))
		 	throw std::invalid_argument("file "+filename+ " not found in GroupClassification::GroupClassification");		
		ifstream is{p};	
		return from_stream(is);
	}
	static GroupClassification from_stream(istream& is) {
		GroupClassification result;
		string line;
		while (std::getline(is,line)) {						
			result.add(make_unique<AbstractLieGroup<true>>(line.c_str(),symbols()));
		}
		return result;
	}
};

GroupClassification nilpotent_lie_groups_3() {
	return GroupClassification::from_disk("classifications/nilpotent3.list");
}
GroupClassification nilpotent_lie_groups_4() {
	return GroupClassification::from_disk("classifications/nilpotent4.list");
}
GroupClassification nilpotent_lie_groups_5() {
	return GroupClassification::from_disk("classifications/nilpotent5.list");
}
GroupClassification nilpotent_lie_groups_6() {
	return GroupClassification::from_disk("classifications/nilpotent6.list");
}
GroupClassification nilpotent_lie_groups_7() {
	return GroupClassification::from_disk("classifications/nilpotent7.list");	
}
GroupClassification nice_nilpotent_lie_groups_6() {
	return GroupClassification::from_disk("classifications/nicenilpotent6.list");
}
GroupClassification nice_nilpotent_lie_groups_7() {
	return GroupClassification::from_disk("classifications/nicenilpotent7.list");
}
GroupClassification nice_nilpotent_lie_groups_8() {
	return GroupClassification::from_disk("classifications/nicenilpotent8.list");
}
GroupClassification nice_nilpotent_lie_groups_9() {
	return GroupClassification::from_disk("classifications/nicenilpotent9.list");
}



/* Nonnice nilpotent Lie groups of dimension 7, as per
 	Diego Conti, Federico Rossi. Construction of nice nilpotent Lie groups. Journal of Algebra, (2019) 525:311-340. doi:10.1016/j.jalgebra.2019.01.020 [Table 2]
 */
class NonniceNilpotentLieGroups7 : public GroupClassification {
public:
	NonniceNilpotentLieGroups7() {
		//reducible
		add("0,0,12,13,0,14+23+25,0");


		add("0,0,12,13,14,15+23,16+23+24");
		add("0,0,12,13,14,15+23,16+24+25-34");		
		add("0,0,12,13,14+23,15+24,16+23+25");	
		add("0,0,12,13,14+23,15+24,-16+23-25");	//123457H1

		add("0,0,12,13,14+23,15+24,23");
		add("0,0,12,13,14+23,25-34,23");
		add("0,0,12,13,14,23,16+25+24-34");
		add("0,0,12,13,14,23,15+25+26-34");
		add("0,0,12,13,23,15+24,14+16+25+34");
		add("0,0,12,13,23,15+24,14+16-25+34");
		add("0,0,12,13,23,15+24,14+16+34");
		add("0,0,12,13,23,15+24,14+16+[lambda]*25+26+34-35",N.lambda);

		add("0,0,12,13,23,-14-25,16+25-35");
		add("0,0,12,13,23,-14-25,15+16+24+[lambda]*25-35",N.lambda);
		add("0,0,12,13,14,0,15+23+26");
		add("0,0,12,13,14+23,0,15+24+26");
		add("0,0,12,13,0,14+25,25+35+16");  	//scritta male in construction
		add("0,0,12,13,0,14+23+25,16+24+35");
		add("0,0,12,13,0,14+23+25,26-34"); 
		add("0,0,12,13,0,14+23+25,15+26-34");
		add("0,0,0,12,14+23,15-34,16+23-35");
		add("0,0,12,13,0,25+23,14");
		add("0,0,12,13,0,14+23,23+25");
		add("0,0,12,0,13,23+24,15+16+25+[lambda]*26+34",N.lambda);
		add("0,0,12,13,0,14+23+25,0");
		add("0,0,12,13,0,14+25+23,15");
		add("0,0,0,12,14+23,23,15-34");
		add("0,0,12,0,23,14,16+26+25-34");
		add("0,0,12,0,24+13,14,15+23+1/2*26+1/2*34");
		add("0,0,12,0,24+13,14,15+[lambda]*23+34+46",N.lambda);
		add("0,0,0,12,13,14+24-35,25+34");
		add("0,0,0,12,13,15+35,25+34");
		add("0,0,0,12,23,-13,15+16+26-2*34");
		add("0,0,0,12,23,-13,[-lambda]*16+[lambda]*25+2*26-2*34",N.lambda);
		add("0,0,0,12,14+23,0,15-34+36");
		add("0,0,0,12,14+23,0,15-34+24+36");
		add("0,0,12,0,0,13+14,15+23");
		add("0,0,12,0,0,2*13+14+25,15+2*23-24");	//cambiate le costanti in modo che a_iR agisca in modo antisimmetrico
	}
};


#endif
