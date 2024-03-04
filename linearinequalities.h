#ifndef LINEARINEQUALITIES_H
#define LINEARINEQUALITIES_H
 
using namespace GiNaC;

template<typename iter_begin,typename iter_end, typename Closure>
iter_begin minimize_value(iter_begin begin, iter_end end, Closure&& closure) {
	int minimum=numeric_limits<int>::max();
	auto minimum_iter=begin;
	while (begin!=end) {
		int value=closure(*begin);
		if (value<minimum) {
			minimum=value;
			minimum_iter=begin;
		}
		++begin;
	}
	return minimum_iter;
}

template<typename iter_begin,typename iter_end, typename T, typename Closure>
T transform_accumulate(iter_begin begin, iter_end end, T value, Closure&& closure) {
	while (begin!=end) value+=closure(*begin++);
	return value;
}

struct ComparableToVariable {
	list<ex> compare_to_x;
	list<ex> x_compares_to;
	ComparableToVariable subs(const lst& subs) const {
		auto result=*this;
		for (auto& x: result.compare_to_x) x=x.subs(subs);
		for (auto& x: result.x_compares_to) x=x.subs(subs);
		return result;
	}	
};

ostream& operator<<(ostream& os, const ComparableToVariable& c) {
	return os<<"ComparableToVariable:  "<<horizontal(c.x_compares_to)<<"<x<"<<horizontal(c.compare_to_x);
}

class OutOfMemory : public std::runtime_error {
public:
	OutOfMemory(string s) : std::runtime_error(s) {}
};



//inequalities 
template<typename Compare, typename NegateCompare>
class GenericLinearInequalities {
	list<ex> inequalities;	//inequality of the form x Compare 0, with the negation being 0 NegateCompare x.
	Compare comparator;
	NegateCompare negative_comparator;
public:
	template<typename Iter1, typename Iter2>
	void insert(Iter1 begin, Iter2 end) {
		inequalities.insert(inequalities.end(),begin,end);
	}
	void insert(ex x) {
		inequalities.push_back(x);
	}
	int count_occurrences(ex variable) const {
		auto present = [variable](ex inequality) {
			exset found;
			return inequality.find(variable,found)? 1 : 0;
		};
		return transform_accumulate(inequalities.begin(),inequalities.end(),0,present);
	}
	bool remove_constant_inequalities() {
		auto iter=inequalities.begin();
		while (iter!=inequalities.end()) {
			if (comparator(iter->evalf(),0)) iter=inequalities.erase(iter);
			else if (negative_comparator(iter->evalf(),0)) return false;
			else ++iter;	//this is the case when *iter depends on a parameter
		}
		return true;
	}

	ComparableToVariable eliminate_variable(ex x) {
		if (inequalities.size()>10000) throw OutOfMemory("too many inequalities: "+to_string(inequalities.size()));		
		ComparableToVariable result;
		for (auto it=inequalities.begin();it!=inequalities.end();) {
			ex coeff=it->coeff(x);
			assert(is_a<numeric>(coeff));
			if (coeff.is_zero()) ++it;
			else if (coeff>0) {				
				result.x_compares_to.push_back(x-*it/coeff);
				it=inequalities.erase(it);
			}			
			else {				
				result.compare_to_x.push_back(x-*it/coeff);
				it=inequalities.erase(it);
			}		
		}		
		return result;
	}
	template<typename Variable, typename Container>
	void GetSymbols(Container& container) const {
		Wedge::GetSymbols<Variable>(container,inequalities.begin(),inequalities.end());
	}
	bool empty() const {
		return inequalities.empty();
	}
	list<ex> list_of_inequalities() const {
		return inequalities;
	}
	void eliminate_duplicates() {
		inequalities.sort(ex_is_less());
		inequalities.unique();
	}
	void subs(const lst& subs) {		
		for (auto& x: inequalities)
		 	x=x.subs(subs);		
		eliminate_duplicates();
	}
};

struct GreaterThan {
	template<typename U,typename T>
	bool operator() (U x, T y) const {return x>y;}
};
struct GreaterOrEqualThan {
	template<typename U,typename T>
	bool operator() (U x, T y) const {return x>=y;}
};
struct LessThan {
	template<typename U,typename T>
	bool operator() (U x,T  y) const {return x<y;}
};
struct LessOrEqualThan {
	template<typename U,typename T>
	bool operator() (U x,T y) const {return x<=y;}
};



struct IntermediateInequalities {
	exvector variables_in_elimination_order;
	vector<ComparableToVariable> strictly_comparable, comparable;
	ex infinity=symbol("inf");
	
	ex ubound(const ComparableToVariable& comparable) const {
		return comparable.compare_to_x.empty()? infinity : *std::min_element(comparable.compare_to_x.begin(),comparable.compare_to_x.end());
	}
	ex lbound(const ComparableToVariable& comparable) const {
		return comparable.x_compares_to.empty()? -infinity : *std::max_element(comparable.x_compares_to.begin(),comparable.x_compares_to.end());
	}
	ex pick_element_less_than(ex ubound) const {
		return ubound==infinity? 0 : ubound-1;
	}
	ex pick_element_in_interval_bounded_below(ex lbound, ex ubound) const {		
		return ubound==infinity? lbound+1 : (ubound+lbound)/2;
	}
	ex pick_element_in_interval(ex lbound, ex ubound) const {		
		if (lbound==-infinity) return pick_element_less_than(ubound);
		else return pick_element_in_interval_bounded_below(lbound, ubound);
	}
	ex find_solution(const ComparableToVariable& strictly_comparable, const ComparableToVariable& comparable) const {	//assumes constant values
		auto strict_ubound=this->ubound(strictly_comparable);
		auto ubound=this->ubound(comparable);
		auto strict_lbound=this->lbound(strictly_comparable);
		auto lbound=this->lbound(comparable);		
//		cout<<lbound<<"<= x <="<<ubound<<endl;
//		cout<<strict_lbound<<"< x <"<<strict_ubound<<endl;
		if (lbound==ubound) return ubound;
		//otherwise, we may assume lbound< ubound, lbound < strict_ubound		
		if (lbound>strict_lbound || strict_lbound==-infinity) strict_lbound=lbound;
		if (ubound<strict_ubound || strict_ubound==infinity) strict_ubound=ubound;
//		cout<<lbound<<"<= x <="<<ubound<<endl;
//		cout<<strict_lbound<<"< x <"<<strict_ubound<<endl;
		return pick_element_in_interval(strict_lbound,strict_ubound);
	}

public:
	void update(ex variable, const ComparableToVariable& strictly_comparable, const ComparableToVariable& comparable) {
		variables_in_elimination_order.push_back(variable);
		this->strictly_comparable.push_back(strictly_comparable);
		this->comparable.push_back(comparable);
	}
	lst solution() const {
		lst result;		
		auto s=strictly_comparable.rbegin();
		auto c=comparable.rbegin();
		for (auto v=variables_in_elimination_order.rbegin();v!=variables_in_elimination_order.rend();++v,++s,++c) {
/*			cout<<"strictly: "<<*s<<endl;
			cout<<"nonstrictly: "<<*c<<endl;			
			cout<<"variable: "<<*v<<endl;
			cout<<"strictly subst: "<<s->subs(result)<<endl;
			cout<<"nonstrictly subst: "<<c->subs(result)<<endl;			

*/
			result.append(*v== find_solution(s->subs(result),c->subs(result)));
//			cout<<"sol="<<result<<endl;
		}
		return result;
	}
	string dbg_print() const {	
		stringstream summary;
		auto strictly=strictly_comparable.begin();
		auto comp=comparable.begin();		
		for (auto x: variables_in_elimination_order) {
			stringstream s;
			s<<"variable "<<x<<": ";
			s<<*strictly++<<endl;
			summary<<">0: size "<<s.str().size();
			s.str();
			s<<*comp++<<endl;
			summary<<", \\geq0: size "<<s.str().size();
		}
		return summary.str();
	}
};


ostream& operator<<(ostream& os, const IntermediateInequalities& o) {
	return os<<"IntermediateInequalities object with sizes : "<<o.dbg_print();
}

template<typename Variable>
class LinearInequalities {
	GenericLinearInequalities<GreaterThan,LessOrEqualThan> positive_;
	GenericLinearInequalities<GreaterOrEqualThan,LessThan> nonnegative_;
	lst subs;
	list<ex> unknowns;	
	IntermediateInequalities intermediate_inequalities;

	ex variable_that_appears_in_the_least_equations() const {
		return *minimize_value(unknowns.begin(),unknowns.end(),[this] (ex variable) {
			int occurrences=positive_.count_occurrences(variable)+nonnegative_.count_occurrences(variable);
			return occurrences? occurrences : numeric_limits<int>::max();
		});
	}	
	
	void eliminate() {
		if (!unknowns.empty())
			eliminate(variable_that_appears_in_the_least_equations());
	}

	template<typename LinearInequalities> 
	void insert_differences(list<ex> bigger, list<ex> smaller, LinearInequalities& destination) {
		bigger.sort(ex_is_less());
		bigger.unique();
		smaller.sort(ex_is_less());
		smaller.unique();		
		for (auto big: bigger)
		for (auto small: smaller)
			destination.insert((big-small).expand());
	}
	void eliminate(ex x) {
		ComparableToVariable strictly_comparable=positive_.eliminate_variable(x);
		ComparableToVariable comparable=nonnegative_.eliminate_variable(x);
		intermediate_inequalities.update(x,strictly_comparable,comparable);
		insert_differences(strictly_comparable.compare_to_x,strictly_comparable.x_compares_to,positive_);
		insert_differences(strictly_comparable.compare_to_x,comparable.x_compares_to,positive_);
		positive_.eliminate_duplicates();
		insert_differences(comparable.compare_to_x,strictly_comparable.x_compares_to,positive_);
		insert_differences(comparable.compare_to_x,comparable.x_compares_to,nonnegative_);
		nonnegative_.eliminate_duplicates();	
	}
	void update_symbols() {
		set<ex,ex_is_less> variables;
		positive_.GetSymbols<Variable>(variables);
		nonnegative_.GetSymbols<Variable>(variables);
		GetSymbols<Variable>(variables,subs);		
		unknowns=list<ex>{variables.begin(),variables.end()};
	}
	bool remove_constant_inequalities() {
		return positive_.remove_constant_inequalities() && nonnegative_.remove_constant_inequalities();
	}
	bool empty() const {
		return positive_.empty() && nonnegative_.empty();
	}
	template<typename Iterator>
	pair<lst,lst> eqns_and_symbols(Iterator begin, Iterator end) const {
		lst eqns=subs;
		list<ex> symbols;
		for (auto i=begin;i!=end;++i)
			eqns.append(i->subs(subs)==0);
		GetSymbols<Variable>(symbols,eqns.begin(),eqns.end());
		return make_pair(eqns,lst{symbols});
	}
	static LinearInequalities impossible_system() {
		ex minus_one=-1;
		return LinearInequalities{}.positive(&minus_one,&minus_one+1);
	}
	void assign_values_to_free_variables(lst& solution) const {
	/*	every entry of subs has one of the following forms

		x=x, where x appears in the inequalities, so it takes a value in the solution
		x=x, where x does not appear in the inequalities, so it should be assigned a value
		x=a_1y_1+...+a_ny_n, where the y_i may appear in the inequalities or not

		this function updates solution by assigning zero to variables of the second type

	*/
		for (auto eq : subs)
			if (eq.lhs()==eq.rhs()) {
				ex equation_with_substituted_parameters=eq.subs(solution);
				if (is_a<Variable>(equation_with_substituted_parameters.lhs()))
					solution.append(equation_with_substituted_parameters.lhs()==0);
			}
	}

	//applies the Fourier-Motzkin algorithm to determine whether the inequalities have a common solution
	bool destructive_has_solution() {
		try {
				if (!remove_constant_inequalities()) return false;
			while (!empty()) {
				if (!remove_constant_inequalities()) return false;
				eliminate();
			}
		}
		catch (const OutOfMemory& p) {
			cerr<<"out of memory when solving linear inequalities (try a better algorithm!) "<<p.what()<<endl;
			return false;
		}
		return true;
	}
	
	lst destructive_find_solution()  {
		if (!destructive_has_solution()) return  lst{};
		auto solution_to_inequalities=intermediate_inequalities.solution() ;
		assign_values_to_free_variables(solution_to_inequalities);
		lst result=solution_to_inequalities;
		for (auto eq : subs)
			if (eq.lhs()!=eq.rhs())
				result.append(eq.subs(solution_to_inequalities));
		return result;
	}

public:
	template<typename Iterator>
	LinearInequalities positive(Iterator begin, Iterator end) const {
		LinearInequalities result=*this;
		result.positive_.insert(begin,end);
		result.positive_.subs(subs);		
		result.update_symbols();
		return result;
	}
	template<typename Iterator>
	LinearInequalities nonnegative(Iterator begin, Iterator end) const {
		LinearInequalities result=*this;
		result.nonnegative_.insert(begin,end);
		result.nonnegative_.subs(subs);		
		result.update_symbols();
		return result;
	}
	template<typename Iterator>
	LinearInequalities zero(Iterator begin, Iterator end) const {
		LinearInequalities result=*this;
		auto eqns_and_symbols=this->eqns_and_symbols(begin,end);
		result.subs =ex_to<lst>(lsolve(eqns_and_symbols.first,eqns_and_symbols.second));
		if (eqns_and_symbols.first.nops() && !result.subs.nops())
			return impossible_system(); //there is no solution to the equalities
		result.positive_.subs(result.subs);
		result.nonnegative_.subs(result.subs);
		result.update_symbols();
		return result;
	}
	LinearInequalities	positive(const list<ex>& l) const {
		return positive(l.begin(),l.end());
	}
	LinearInequalities	nonnegative(const list<ex>& l) const {
		return nonnegative(l.begin(),l.end());
	}
	LinearInequalities zero(const list<ex>& l) const {
		return zero(l.begin(),l.end());
	}

	bool has_solution() const {
		auto copy=*this;
		return move(copy).destructive_has_solution();
	}
	lst find_solution() const {
		auto copy=*this;
		return move(copy).destructive_find_solution();
	}

	list<ex> get_positive() const {return positive_.list_of_inequalities();}
	list<ex> get_nonnegative() const {return nonnegative_.list_of_inequalities();}


	list<ex> get_zero() const {
		list<ex> result;
		for (auto x: subs)
			if (x.lhs()!=x.rhs()) result.push_back(x.lhs()-x.rhs());
		return result;
	}
};

template<typename T>
ostream& operator<<(ostream& os, const LinearInequalities<T>& l) {
	auto positive=l.get_positive();
	auto nonnegative=l.get_nonnegative();
	auto zero=l.get_zero();
	if (positive.empty() && nonnegative.empty() && zero.empty()) return os<<"true";
	string sep;
	if (!positive.empty()) {os<<horizontal(positive)<<">0"; sep=", ";}
	if (!nonnegative.empty()) {os<<sep<<horizontal(nonnegative)<<"\\geq0"; sep=", ";}
	if (!zero.empty()) os<<sep<<horizontal(zero)<<"=0";
	return os;
}

//represent a list of linear inequalities logically separated by an or, e.g. x<0 || x>0
template<typename Variable>
class AlternativeLinearInequalities {
	list<LinearInequalities<Variable>> alternatives;
public:	
	AlternativeLinearInequalities operator||(const AlternativeLinearInequalities& ineq) const {
		auto result=*this;
		result.alternatives.insert(result.alternatives.end(),ineq.alternatives.begin(),ineq.alternatives.end());
		return result;
	}
	AlternativeLinearInequalities operator||(const LinearInequalities<Variable>& ineq) const {
		auto result=*this;
		result.alternatives.emplace_back(ineq);
		return result;
	}
	AlternativeLinearInequalities operator&&(const AlternativeLinearInequalities<Variable>& other) const {
		AlternativeLinearInequalities result;
		for (auto& ineq : alternatives)
		for (auto& ineq2 : other.alternatives)
			result.alternatives.push_back(
				ineq.positive(ineq2.get_positive()).nonnegative(ineq2.get_nonnegative()).zero(ineq2.get_zero())
			);
		return result;
	}
	AlternativeLinearInequalities operator&&(const LinearInequalities<Variable>& ineq) const {
		AlternativeLinearInequalities result;
		for (auto& x: alternatives)
			result.alternatives.push_back(x.positive(ineq.get_positive()).nonnegative(ineq.get_nonnegative()).zero(ineq.get_zero()));
		return result;
	}
	static AlternativeLinearInequalities positive(ex x)  {
		AlternativeLinearInequalities result;	
		result.alternatives.push_back(LinearInequalities<Variable>{}.positive(&x,&x+1));
		return result;
	}
	static AlternativeLinearInequalities negative(ex x)  {
		return positive(-x);		
	}
	static AlternativeLinearInequalities nonnegative(ex x)  {
		AlternativeLinearInequalities result;
		result.alternatives.push_back(LinearInequalities<Variable>{}.nonnegative(&x,&x+1));
		return result;
	}
	static AlternativeLinearInequalities nonpositive(ex x)  {
		return nonnegative(-x);
	}
	static AlternativeLinearInequalities nonzero(ex x)  {
		return positive(x) || negative(x);
	}
	static AlternativeLinearInequalities zero(ex x)  {
		AlternativeLinearInequalities result;
		result.alternatives.push_back(LinearInequalities<Variable>{}.zero(&x,&x+1));
		return result;
	}
	bool has_solution() const {
		return any_of(alternatives.begin(), alternatives.end(),
			[] (auto& ineq) {return ineq.has_solution();});
	}
	lst find_solution() const {
		for (auto& ineq : alternatives) {			
			auto sol=ineq.find_solution();
			if (sol.nops()) return sol;
		}
		return {};
	}
	auto alternative_begin() const {return alternatives.begin();}
	auto alternative_end() const {return alternatives.end();}
};

template<typename Variable>
ostream& operator<<(ostream& os, AlternativeLinearInequalities <Variable>& alts) {
	bool first=true;
	for (auto i=alts.alternative_begin();i!=alts.alternative_end();++i) {
		if (!first) os<<"\t|| ";
		os<<*i<<endl;
		first=false;		
	}
	return os;
}

namespace inequalities {
	
AlternativeLinearInequalities<symbol> positive(ex x) {return AlternativeLinearInequalities<symbol>::positive(x);}
AlternativeLinearInequalities<symbol> nonpositive(ex x) {return AlternativeLinearInequalities<symbol>::nonpositive(x);}
AlternativeLinearInequalities<symbol> negative(ex x) {return AlternativeLinearInequalities<symbol>::negative(x);}
AlternativeLinearInequalities<symbol> nonnegative(ex x) {return AlternativeLinearInequalities<symbol>::nonnegative(x);}
AlternativeLinearInequalities<symbol> zero(ex x) {return AlternativeLinearInequalities<symbol>::zero(x);}
AlternativeLinearInequalities<symbol> nonzero(ex x) {return AlternativeLinearInequalities<symbol>::nonzero(x);}

}

void test_linear_inequalities() {
	symbol x("x"),y("y"),z("z");
	exvector positive={x+y,x+z};
	exvector zero={z+x+y};
	exvector nonnegative={};
	auto l=LinearInequalities<symbol>{}.positive(positive.begin(),positive.end()).nonnegative(nonnegative.begin(),nonnegative.end()).zero(zero.begin(),zero.end());
	cout<<l<<endl;
	cout<<l.find_solution()<<endl;	
}

void test_linear_inequalities2() {
	symbol x("x"),y("y"),z("z");
	exvector positive={x+y,-2*x+z};
	exvector zero={z,3*x+2*y};
	exvector nonnegative={y};
	auto l=LinearInequalities<symbol>{}.positive(positive.begin(),positive.end()).nonnegative(nonnegative.begin(),nonnegative.end()).zero(zero.begin(),zero.end());
	cout<<l<<endl;
	cout<<l.find_solution()<<endl;	
}

void test_alternative_linear_inequalities() {
	using namespace inequalities;
	symbol x("x"),y("y"),z("z");
	auto p=positive(x+y)&& zero(z+x+y) && nonnegative(y) && zero(y+1);
	auto q=p || negative(y);
	cout<<q<<endl;
	cout<<q.find_solution()<<endl;

}

void test_double_free() {
	symbol w1("w1"),w2("w2"),w3("w3"),w4("w4"),w5("w5"),w6("w6"),w7("w7"),w8("w8"),w9("w9");
	list<ex> positive={-w9+w3+w7,2* w5-w8,w2,w3+w7-w8,w2-w9+w8,2 *w5-w9,w6+w4-w8,w1,w9+w1-w8,-w9+w6+w4};
	list<ex> nonnegative={ -w2+w9-w3,-w5-w2+w7,-w5+w9-w4,-w7+w8,-w3-w4+w8,-w2-w6+w8,-w6+w7,w6-w3-w1,w5-w4,-w2+w4-w1,w9-w7-w1,w5-w4-w1,-w3+w4,-w5-w1+w8,-w2+w3,w2-w1,w9-w8,-w2+w6-w4,-w5+w6,-w6+w7-w1};
	auto p=LinearInequalities<symbol>{}.positive(positive).nonnegative(nonnegative);
	cout<<p.has_solution()<<endl;
}

#endif
