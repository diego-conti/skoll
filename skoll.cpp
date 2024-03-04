#include <ratatoskr/ratatoskr.h>
#include <numeric>
#include <optional>
#include "linearsolve.h"
#include "derivations.h"
#include "classification.h"
#include "torusinder.h"
#include "graded.h"
#include "linearextensions.h"
#include "linearinequalities.h"
#include "filtered.h"
#include "antidiagonal.h"
#include "sigmadiagonal.h"


matrix ricci_tensor(const Manifold& G, matrix metric_on_frame) {
	auto P=PseudoRiemannianStructureByMatrix::FromMatrixOnFrame (&G,G.e(),metric_on_frame);
	PseudoLeviCivitaConnection omega(&G,P);
	matrix ricci_tensor=omega.RicciAsMatrix();
	assert(ricci_tensor.cols()==G.Dimension());
	assert(ricci_tensor.rows()==G.Dimension());
	return ricci_tensor;
}
   
struct MetricAndRicci {	
	matrix g;
	ex ric;
	ex Ric;	
	int score=std::numeric_limits<int>::max();
	MetricAndRicci()=default;	
	MetricAndRicci(const LieGroup& G, const matrix& g) : g{g}, ric {ricci_tensor(G,g)}, Ric{(g.inverse()*ric).evalm()} {
		static CocoaPolyAlgorithms::Initializer initialize_cocoa;
		set<ex,ex_is_less> entries(begin(Ric),end(Ric));				
		GetSymbols<symbol>(symbols,entries.begin(),entries.end());
		for (auto& x: entries) if (!x.is_zero()) ideal.push_back(x.numer());
		score=ideal.size();
		try {
		if (CocoaPolyAlgorithms_R::RadicalContains<symbol>(symbols, ideal.begin(),ideal.end(),g.determinant())) score+=100;
		}
		 catch (const CoCoA::ErrorInfo& err)
		 {
			cerr<<err<<endl;
			throw err;
		 }  
	}
	exvector reduced_ideal() const {		
		return CocoaPolyAlgorithms_R::IdealReduce<symbol>(symbols, ideal.begin(),ideal.end());
	}
	bool no_solution() const {return score>=100;}
private:
	exvector symbols;
	exvector ideal;
};


string to_pair(int i,int j, int bound) {
	if (bound<10) return to_string(i)+to_string(j);
	else return to_string(i)+","+to_string(j);
}

//assumes g is symmetric
string nonzero_entries(const matrix& g) {
	stringstream s;
	for (int i=0;i<g.rows();++i)
	for (int j=i+1;j<g.rows();++j)
		if (!g(i,j).expand().is_zero()) s<<to_pair(i+1,j+1,g.rows())<<"\\ ";
	if (s.str().empty()) return "-";	//empty list
	return s.str();
}

ostream& operator<<(ostream& os, const MetricAndRicci& m ) {
	if (m.no_solution()) return os;
	os<<nonzero_entries(m.g)<<endl;
	if (m.score) 
		os<<"if "<<horizontal(m.reduced_ideal())<<"=0 "<<m.g<<endl;
	return os;
}


struct FindMetricParameters {
	vector<matrix> imaginary_derivations_in_torus;
	optional<exvector> grading;
	FindMetricParameters(const exvector& v) : grading{v} {}
	FindMetricParameters(const TorusInDer& t) : imaginary_derivations_in_torus{t.imaginary_derivations_in_torus()} , grading{t.grading()} {}	
};

MetricAndRicci best_sigmadiagonal_metric(const LieGroup& G, const FindMetricParameters&) {
	MetricAndRicci best;		
	for (auto g: SigmaDiagonalMetrics{G.Dimension()}) {			
		MetricAndRicci metric_and_ricci(G,g);
		if (metric_and_ricci.score==0) 
			return metric_and_ricci;				
		else if (metric_and_ricci.score<best.score) best=metric_and_ricci;
	}	
	return best;
}
bool find_sigmadiagonal_metric(const LieGroup& G, const FindMetricParameters& p,ostream& os) {
	auto best=best_sigmadiagonal_metric(G,p);
	os<<best<<"\\\\"<<endl;
	return !best.no_solution();
}


bool find_filtered_metric(const LieGroup& G,  const FindMetricParameters& p, ostream& os) {	
	for (Filtration f{G}; f; ++f) {
		os<<horizontal(f.basis())<<"&"<<horizontal(f.weights());
		os<<"\\\\"<<endl;			
		return true;
	}	
	os<<"\\\\"<<endl;			
	return false;
}

exvector adapted_basis_from_indices(const exvector& e, const vector<int>& indices) {
	exvector result;
	transform(indices.begin(),indices.end(),back_inserter(result),[e] (int i) {return e[i];});
	return result;
}

bool find_foad_metric(const LieGroup& G, const FindMetricParameters& p,ostream& os) {	
	if (p.grading) {	
		auto sequences=WeightSequencesRespectingOrder{p.grading.value()};	
		for (auto H : sequences) {
			os<<horizontal(H)<<"&"<<horizontal(adapted_basis_from_indices(G.e(),sequences.as_indices(H)));
			os<<"\\\\"<<endl;			
			return true;		
		}	
		os<<"\\text{weights of split torus: } "<<horizontal(p.grading.value());
	}
	else os<<"\\text{ERROR: split torus acts nondiagonally}";
	os<<"\\\\"<<endl;			
	return false;
}

bool find_any_ricciflat_metric(const LieGroup& G, const FindMetricParameters& p,ostream& os) {	
	stringstream foad_stream, filtered_stream, sigma_stream;
	if (find_foad_metric(G,p,foad_stream)) {
		os<<"(G1)--(G5):"<<foad_stream.str();			
		return true;
	}
	else if (find_filtered_metric(G,p,filtered_stream)) {
			os<<"(F1)--(F5):"<<filtered_stream.str();
			return true;
	}
	else if (find_sigmadiagonal_metric(G,p,sigma_stream)) {
		os<<"&"<<sigma_stream.str();
		return true;
	}
	else {
		os<<"no Ricci-flat metric found!\\\\"<<endl;
		return false;
	}
}

string canonical_print_no_brackets(const LieGroup& G) {
	stringstream s;
	s<<latex;
	G.canonical_print(s);
	auto as_string=s.str();
	assert(as_string.front()=='(' && as_string.back()==')');
	return as_string.substr(1,as_string.size()-2);
}

template<typename FindFunction, typename Filter>
void print_table_row(const LieGroup& G,ostream& os, FindFunction& find_metric,int columns_for_lie_algebra, Filter filter) {	
	TorusInDer a{G};	
	if (filter(a))	{
		if (columns_for_lie_algebra>1)		os<<"\\multicolumn{"<<columns_for_lie_algebra<<"}{L}{";
		os<<canonical_print_no_brackets(G);		
		if (columns_for_lie_algebra>1) os<<"}"<<"\\\\ ";		
		os<<"&";
		if (!a.is_a_direct_sum()) os<<"WARNING: torus is not a direct sum of symmetric and skew-symmetric matrices; ";
		//a.print_table_row(os);			
		find_metric(G,a,os);		
	}
}
template<typename FindFunction>
void print_table_row(const LieGroup& G,ostream& os, FindFunction& find_metric,int columns_for_lie_algebra) {	
	auto no_filter=[](auto& ) {return true;};
	print_table_row(G,os,find_metric,  columns_for_lie_algebra, no_filter);
}

void choose_basis_if_one_dimensional(ExVector& v) {
	list<ex> symbols;
	GetSymbols<symbol>(symbols,v.begin(),v.end());
	if (symbols.size()==1) 
		for (auto& x: v)
			x=x.subs(symbols.front()==1);
}

template<typename FindFunction, typename Filter>
void print_table_row_nice(const LieGroup& G,ostream& os, FindFunction& find_metric, int columns_for_lie_algebra,Filter filter) {	
	GL gl(G.Dimension());	
	if (columns_for_lie_algebra>1)		os<<"\\multicolumn{"<<columns_for_lie_algebra<<"}{L}{";
	os<<canonical_print_no_brackets(G);		
	if (columns_for_lie_algebra>1) os<<"}"<<"\\\\ ";		
	os<<"&";
	auto der=diagonal_derivations_on_nice_lie_algebra(G);
	choose_basis_if_one_dimensional(der);
	find_metric(G,der,os);	
}
template<typename FindFunction>
void print_table_row_nice(const LieGroup& G,ostream& os, FindFunction& find_metric,int columns_for_lie_algebra) {	
	auto no_filter=[](auto& ) {return true;};
	print_table_row_nice(G,os,find_metric, columns_for_lie_algebra, no_filter);	
}

enum class ClassOfLieAlgebras {
	NICE, NONNICE, ALL
};

struct Parameters {
	unique_ptr<LieGroup> G;
	int d;
	ratatoskr::GlobalSymbols symbols;
	ClassOfLieAlgebras class_of_lie_algebras=ClassOfLieAlgebras::ALL;
	int columns_for_lie_algebra=1;
};

auto parameter_description= ratatoskr::make_parameter_description(
		ratatoskr::alternative("lie algebra or dimension of lie algebras to study")(
			"lie-algebra", "lie algebra, possibly with parameters", ratatoskr::lie_algebra(&Parameters::G, &Parameters::symbols)
		)(
			"dimension", "study all (nice) Lie algebras of dimension (3 through 9 are allowed)", &Parameters::d
		),
		ratatoskr::alternative("nice|non-nice|all")(
			"non-nice", "only study non-nice Lie algebras (implemented for dimension 7)",ratatoskr::generic_option(&Parameters::class_of_lie_algebras, [] () {return ClassOfLieAlgebras::NONNICE;})
		)(
			"nice", "only study nice Lie algebras (using diagonal torus)",ratatoskr::generic_option(&Parameters::class_of_lie_algebras, [] () {return ClassOfLieAlgebras::NICE;})
		)(
					
			"all", "all Lie algebras",ratatoskr::generic_option(&Parameters::class_of_lie_algebras, [] () {return ClassOfLieAlgebras::ALL;})
		),
		"columns","columns to use to represent the Lie algebra in the output when printing a table",&Parameters::columns_for_lie_algebra
	);

template<typename FindFunction>
void study_one(const LieGroup& G, ostream& os,  FindFunction& f) {
	print_table_row(G,os,f,1);
}

template<typename Classification, typename... FindFunctionAndFilter>	
void study_classification(Parameters& parameters, ostream& os, const Classification& classification, FindFunctionAndFilter... f) {
	if (parameters.class_of_lie_algebras==ClassOfLieAlgebras::NICE) 
		for (auto& G: classification)
			print_table_row_nice(*G,os,f...);	
	else	
		for (auto& G: classification)
			print_table_row(*G,os,f...);		
}

template<typename... FindFunctionAndFilter>
void study_all(Parameters& parameters, ostream& os, FindFunctionAndFilter... f) {
		os<<"%\\begin{array}{ccc}"<<endl;		
		if (parameters.class_of_lie_algebras==ClassOfLieAlgebras::NONNICE) {
			if (parameters.d!=7) 
				cerr<<"not implemented"<<endl;
			else study_classification(parameters,os, NonniceNilpotentLieGroups7(),f...);			
		}
		else 
		switch (parameters.d) {
		case 3:			
			study_classification(parameters, os, nilpotent_lie_groups_3(),f...);
			break;
		case 4:			
			study_classification(parameters, os, nilpotent_lie_groups_4(),f...);
			break;
		case 5:			
			study_classification(parameters, os, nilpotent_lie_groups_5(),f...);
			break;
		case 6:			
			study_classification(parameters, os, 
				parameters.class_of_lie_algebras==ClassOfLieAlgebras::NICE ? nice_nilpotent_lie_groups_6() : nilpotent_lie_groups_6(),
				f...);
			break;
		case 7: {
			study_classification(parameters, os, 
				parameters.class_of_lie_algebras==ClassOfLieAlgebras::NICE ? nice_nilpotent_lie_groups_7() : nilpotent_lie_groups_7(),
				f...);
			break;
		}
		case 8:
			study_classification(parameters, os, nice_nilpotent_lie_groups_8(),f...);			
			break;	
		case 9:			
			study_classification(parameters, os, nice_nilpotent_lie_groups_9(),f...);		
			break;				
		default:
			cerr<<"unsupported dimension "<<parameters.d<<endl;
		}
		os<<"%\\end{array}"<<endl;
}

auto program1=ratatoskr::make_program_description(
	"sigma-diagonal", "study sigma-diagonal metrics", parameter_description, [] (Parameters& parameters, ostream& os) {
		if (parameters.G)
			study_one(*parameters.G,os,find_sigmadiagonal_metric);
		else			
			study_all(parameters,os,find_sigmadiagonal_metric,parameters.columns_for_lie_algebra);		
	}
);	



auto program2=ratatoskr::make_program_description(
	"graded", "study gradings satisfying (G1)--(G5)", parameter_description, [] (Parameters& parameters, ostream& os) {
		if (parameters.G)
			study_one(*parameters.G,os,find_foad_metric);
		else			
			study_all(parameters,os,find_foad_metric,parameters.columns_for_lie_algebra);		
	}
);	

auto program3=ratatoskr::make_program_description(
	"filtered", "study filtrations satisfying (F1)--(F5)", parameter_description, [] (Parameters& parameters, ostream& os) {
		if (parameters.G)
			study_one(*parameters.G,os,find_filtered_metric);
		else			
			study_all(parameters,os,find_filtered_metric,parameters.columns_for_lie_algebra);		
	}
);	

auto program4=ratatoskr::make_program_description(
	"any", "find a Ricci-flat metric of any type", parameter_description, [] (Parameters& parameters, ostream& os) {
		if (parameters.G)
			study_one(*parameters.G,os,find_any_ricciflat_metric);
		else			
			study_all(parameters,os,find_any_ricciflat_metric,parameters.columns_for_lie_algebra);		
	}
);	

void print_derivations(const LieGroup& G, ostream& os) {
	TorusInDer torus{G};
	torus.print(os);
	os<<endl;	
	Grading grading{torus.linear_group(),torus.grading().value(),torus.nilradical()};
	grading.print(os);
}

auto program5=ratatoskr::make_program_description(
	"derivations", "print the radical in the Lie algebra of derivations, decomposed as the sum of a maximal torus and the graded algebra of nilpotent derivations", parameter_description, [] (Parameters& parameters, ostream& os) {
		if (parameters.G) print_derivations(*parameters.G,os);
	}
);


int main(int argv, char** argc) {				
	ratatoskr::alternative_program_descriptions(program1,program2,program3,program4,program5).run(argv,argc);
}

