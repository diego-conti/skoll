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

#include <wedge/wedge.h>
#include <optional>

using namespace GiNaC;
using namespace std;
using namespace Wedge;

ex Xbracket(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A, ex X, ex Y) {
	ex Ax=V.Action<VectorField>(A,X);
	ex Ay=V.Action<VectorField>(A,Y);
	ex Axy=V.Action<VectorField>(A,G.LieBracket(X,Y));
	return G.LieBracket(Ax,Y)+G.LieBracket(X,Ay)-Axy;
}

exvector Xbrackets(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A) {
		exvector Xbrackets;
		for (int i=1;i<=G.Dimension();++i)
		for (int j=i+1;j<=G.Dimension();++j) 
			Xbrackets.push_back(Xbracket(G,V,A,G.e(i),G.e(j)).expand());				
		return Xbrackets;
}

/** Represents a vector space which is sandwiched between a smaller and a larger subspace.
*/
struct VectorSpaceBetween {
	exvector basis_of_smaller_space;
	exvector basis_of_larger_space;
};

/** For a Lie group with parameters, return a VectorSpaceBetween object representing the derivations
	@param G a Lie group of dimension n, with or without parameters
	@param Gl The Lie algebra of GL(n,R), acting on the Lie algebra of g through the identification g=R^n given by the standard coframe of g
	@result A VectorSpaceBetween representing the subspace of Gl corresponding to the space of derivations
	
	The exact space of derivations corresponds to solutions of a linear system depending on parameters. This function computes the space of solutions of a subset of the equations that do not depend on a parameter and the space of elements that satisfy the equations for all values of the parameters
*/	

template<typename Parameter>
VectorSpaceBetween derivations_parametric(const LieGroup& G,const GL& Gl)  {
		auto gl=Gl.pForms(1);
		auto generic_matrix =gl.GenericElement();
		auto X=Xbrackets(G,GLRepresentation<VectorField>(&Gl,G.e()),generic_matrix);
		lst eqns,sol;
		GetCoefficients<VectorField>(eqns,X);
		Wedge::linear_impl::LinearEquationsWithParameters<VectorSpace<DifferentialForm>::Coordinate,Parameter> linear_eqns{eqns,lst{gl.coordinate_begin(),gl.coordinate_end()}};
		linear_eqns.eliminate_linear_equations();
		VectorSpaceBetween result;
		gl.GetSolutionsFromGenericSolution(result.basis_of_larger_space,linear_eqns.solution());
		gl.GetSolutionsFromGenericSolution(result.basis_of_smaller_space,linear_eqns.always_solution());
		return result;
}

//a diagonal matrix (d)_ij is a derivation if d[e_i,e_j]=[de_i,e_j]+[e_i,de_j] i.e. 

template<typename T> ExVector vector_of_symbols(int size, const Name& n) {
	ExVector result;
	result.reserve(size);
	for (int i=1;i<=size;++i) result.emplace_back(T{n(i)});
	return result;
}

//given a nice Lie group G and an index k, return the sequence of pairs (i,j) such that [e_i,e_j] is a nonzero multiple of e_k or nullopt
vector<pair<int,int>> nodes_going_to(const LieGroup& G, OneBased k) {
	vector<pair<int,int>> result;
	ex de_k=G.d(G.e(k));
	if (de_k.is_zero()) return result;
	int n=G.Dimension();
	for (int i=1;i<=n;++i) {
		ex e_i_hook_de_k=Hook(G.e(i),de_k);
		if (!e_i_hook_de_k.is_zero())
			for (int j=i+1;j<=n;++j) {
				if (!TrivialPairing<DifferentialForm>(G.e(j),e_i_hook_de_k).is_zero()) {
					result.push_back(make_pair(i,j));
					break;
				}
			}
	}	
	return result;	
}

//return the equation D(i)+D(j)==D(k)
ex equation(const ExVector& D, int i, int j, int k) {
	return D(i)+D(j)==D(k);
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(DerivationParameter, realsymbol)
ExVector diagonal_derivations_on_nice_lie_algebra(const LieGroup& G) {	
	int n=G.Dimension();
	auto diagonal_derivation=vector_of_symbols<DerivationParameter>(n,N.lambda);
	lst eqns;
	for (int k=1;k<=n;++k) 
		for (auto ij :nodes_going_to(G,k))
		 	eqns.append(equation(diagonal_derivation,ij.first, ij.second, k));	
	auto sol=lsolve(eqns,lst{diagonal_derivation.begin(),diagonal_derivation.end()});
	for (ex& x: diagonal_derivation ) x=x.subs(sol);
	return diagonal_derivation;	
}
