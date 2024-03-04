#include "sumintersect.h"

class TraceScalarProduct {
	const GL& gl;
public:
	TraceScalarProduct(const  GL& gl) : gl{gl} {}
	ex tr(ex gl_element) const {
		return ex_to<matrix>(gl.glToMatrix(gl_element)).trace();
	}
	ex tr_product(ex gl_element1, ex gl_element2) const {
		return ex_to<matrix>((gl.glToMatrix(gl_element1)*gl.glToMatrix(gl_element2)).evalm()).trace();
	}
	VectorSpace<DifferentialForm> perp(const VSpace<DifferentialForm>& subspace, const VSpace<DifferentialForm>& space) const {
		auto gen_elem=space.GenericElement();
		lst eqns;
		for (auto x: subspace.e())
			GetCoefficients<DifferentialForm>(eqns,tr_product(x,gen_elem));
		return space.SubspaceFromEquations(eqns.begin(),eqns.end());
	}
};

VectorSpace<DifferentialForm> nilpotent_derivations(const LieGroup& G, const GL& gl) {
	auto all_der=derivations(G,gl);	
	TraceScalarProduct t{gl};
	auto nil_der=t.perp(all_der,all_der);
	return nil_der;
}


VectorSpace<DifferentialForm> derived_subalgebra(const VectorSpace<DifferentialForm>& subspace, const LieGroup& G) {	 
	set<ex,ex_is_less> generators;
	for (int i=0;i<subspace.Dimension();++i)
	for (int j=i+1;j<subspace.Dimension();++j)	
		generators.insert(G.LieBracket(subspace.e()[i],subspace.e()[j]));
	return VectorSpace<DifferentialForm>{generators.begin(),generators.end()};
}


VectorSpace<DifferentialForm> nilradical_of_linear_algebra(const VectorSpace<DifferentialForm>& linear_algebra, const GL& gl) {	
	auto nilradical=TraceScalarProduct{gl}.perp(linear_algebra,linear_algebra);
    return nilradical;
}

VectorSpace<DifferentialForm> radical_of_linear_algebra(const VectorSpace<DifferentialForm>& linear_algebra, const GL& gl) {
	auto derived_algebra=derived_subalgebra(linear_algebra,gl);
	auto radical=TraceScalarProduct{gl}.perp(derived_algebra,linear_algebra);
    return radical;
}

optional<exvector> diagonal_elements_if_diagonal_matrix(const matrix& m) {
	if (m.rows()!=m.cols()) return nullopt;
	exvector d(m.rows());
	for (int i=0; i<m.rows();++i)
	for (int j=0; j<m.cols();++j) 
		if (!m(i,j).is_zero()) {
			if (i==j) d[i]=m(i,j);
			else return nullopt;
		}
	return d;
}

class SymmetricAndSkewDecomposition {
    const GL& gl;
    VectorSpace<DifferentialForm> space,sym,skew;
    ex transpose(ex x) const {
        matrix m=gl.glToMatrix(x);
        matrix tm=m.transpose();
        return gl.MatrixTo_gl(tm);
    }
public:
    SymmetricAndSkewDecomposition(const GL& gl, const VectorSpace<DifferentialForm>& space) : gl{gl}, space{space} {
        auto x=space.GenericElement();
        auto tx=transpose(x);
        lst skew_eqns,sym_eqns;
        GetCoefficients<DifferentialForm>(skew_eqns,x+tx);
        GetCoefficients<DifferentialForm>(sym_eqns,x-tx);
        skew=space.SubspaceFromEquations(skew_eqns.begin(),skew_eqns.end());
        sym=space.SubspaceFromEquations(sym_eqns.begin(),sym_eqns.end());
    }
    bool is_direct_sum() const {
        return space.Dimension()==sym.Dimension()+skew.Dimension();
    }
    const VectorSpace<DifferentialForm>& symmetric() const {return sym;}
    const VectorSpace<DifferentialForm>& skew_symmetric() const {return skew;}
	const VectorSpace<DifferentialForm>& all() const {return space;}
};

class TorusInDer {
	friend class GradedDerivations;
	GL gl;
	VectorSpace<DifferentialForm> der,n;
    SymmetricAndSkewDecomposition a;
	optional<exvector> weights; //real weights, i.e. weight decomposition using real torus
    exvector imaginary_torus; //basis of the compact torus
	VectorSpace<DifferentialForm> torus() const {
		if (n.Dimension()==der.Dimension()) return {};	//characteristically nilpotent
		auto r=radical_of_linear_algebra(der,gl);    	
    	Subspace<DifferentialForm> n_in_r=r.Subspace(n.e().begin(),n.e().end());
		return n_in_r.Complement();    	
	}
    void print_if_nonempty_or(ostream& os, const exvector& l, string label) const {
        if (l.empty()) os<<label;
        else os<<horizontal(l);
    }
public:
	const GL& linear_group() const {return gl;}
	VectorSpace<DifferentialForm> nilradical() const {return n;}
	optional<exvector> grading() const {return weights;}

	TorusInDer(const LieGroup& G) : gl(G.Dimension()),der{derivations_parametric<LieAlgebraParameter>(G,gl).basis_of_smaller_space}, n{nilradical_of_linear_algebra(der,gl)}, a{SymmetricAndSkewDecomposition{gl,torus()}} {
        auto& aR=a.symmetric();
		auto X=(aR.Dimension()==1? aR.e(1) : aR.GenericElement());
		auto H=gl.glToMatrix(X);
		weights=diagonal_elements_if_diagonal_matrix(H);	
	}

	bool is_a_direct_sum() const {
		return a.is_direct_sum();
	}

	void print(ostream& os) const {
    	os<<"der ="<<horizontal(der.e())<<endl;		
		if (n.Dimension()==der.Dimension()) os<<"characteristically nilpotent"<<endl;
		else {
			os<<"n="<<horizontal(n.e())<<endl;
    		os<<"a_R="<<horizontal(a.symmetric().e())<<endl;			
    		os<<"a_{iR}="<<horizontal(a.skew_symmetric().e())<<endl;
            if (!a.is_direct_sum()) os<<"WARNING: a is not a sum of symmetric and skew-symmetric matrices,"<<horizontal(a.all().e())<<endl;
		}
		if (weights) os<<"weights: "<<horizontal(weights.value());
		else os<<"canonical basis does not diagonalize a_R"<<endl;
	}
	void print_table_row(ostream& os) const {
		if (n.Dimension()==der.Dimension()) 
			os<<"characteristically nilpotent";		
		else {    
			print_if_nonempty_or(os, a.skew_symmetric().e(),"");
            if (!a.is_direct_sum()) os<<" (not a sum, "<<horizontal(a.all().e())<<")";
		}		
		if (weights) os<<"&"<<horizontal(weights.value());
		else os<<"&nondiagonal base; symmetric: "<<horizontal(a.symmetric().e())<<endl;				
	}
	void print_relations(ostream& os) const {
		
	}
	vector<std::initializer_list<OneBased>> relations() const {
		auto& weights=this->weights.value();
		vector<std::initializer_list<OneBased>> result;
		for (int i=0;i<gl.n();++i)
		for (int j=i+1;j<gl.n();++j)
		for (int k=0;k<gl.n();++k)
			if (weights[i]+weights[j]==weights[k]) result.push_back({k+1,i+1,j+1});
		return result;
	}


    vector<matrix> imaginary_derivations_in_torus() const {
        vector<matrix> result;
        auto convert_to_matrix = [this] (ex x) {return gl.glToMatrix(x);};
        transform(a.skew_symmetric().e().begin(),a.skew_symmetric().e().end(),back_inserter(result), convert_to_matrix);
        return result;
    }

	
};



class Grading {
	map<ex,VectorSpace<DifferentialForm>,ex_is_less> g_i;
public:
	Grading(const GL& gl,const exvector& weights) {
		for (int i=0;i<gl.n();++i)
		for (int j=0;j<gl.n();++j)
			g_i[weights[i]-weights[j]].AddGenerator(gl.A(i+1,j+1));
	}
	Grading(const GL& gl,const exvector& weights,const VectorSpace<DifferentialForm>& subspace) : Grading{gl,weights} {		
		map<ex,VectorSpace<DifferentialForm>,ex_is_less> restricted;
		for (auto i_and_g_i : g_i) {
			auto i=i_and_g_i.first;
			auto& g_i=i_and_g_i.second;
			auto h_i=intersected_with(subspace,g_i);
			if (h_i.Dimension()) restricted[i]=h_i;
		}
		swap(restricted,g_i);
	}
	void print(ostream& os) const {
		os<<"Graded space"<<endl;
		for (auto i_and_g_i : g_i) {
			auto i=i_and_g_i.first;
			auto& g_i=i_and_g_i.second;
			os<<i<<": " <<horizontal(g_i.e())<<endl;
		}
	}
};


class GradedDerivations {
	TorusInDer torus;
	Grading grading;
public:
	GradedDerivations(const LieGroup& G) : torus(G), grading{torus.gl,torus.weights.value(),torus.n} {
		grading.print(cout);	
	}
};