WEDGE_DECLARE_NAMED_ALGEBRAIC(WParameter,realsymbol)
 
 //a basis obtained from G.e() by reordering, in such a way that e_i\hook de^j\neq 0 implies i<j and [e(i_hat),e(i)] is in the span of e_n
class OrderedBasis {
    LinearExtension<int> ordered_basis;
    const LieGroup& G;
    static LinearExtension<int> first_compatible_order(const LieGroup& G) {
        vector<int> indices; 
        set<pair<int,int>> poset;
        indices.resize(G.Dimension());
        iota(indices.begin(),indices.end(),0);
        for (int i=0;i<G.Dimension();++i)
        for (int j=0;j<G.Dimension();++j)
            if (!Hook(G.e()[i],G.d(G.e()[j])).is_zero()) poset.emplace(i,j);
        return LinearExtension<int>::begin(indices,poset);
    }
    bool bracket_in_span_of_last(int i, int j) const {
        auto sigma=*ordered_basis;
        ex e_sigma_i=G.e()[sigma[i]];
        ex e_sigma_j=G.e()[sigma[j]];        
        ex e_n=G.e()[sigma.back()];
        return (G.LieBracket(e_sigma_i,e_sigma_j).subs(e_n==0)).expand().is_zero();
    }
    bool is_valid() const {
        for (int i=0;i<G.Dimension();++i)
            if (!bracket_in_span_of_last(i,hat(i))) return false;
        return true;
     }
    void advance_until_valid() {
        while (ordered_basis && !is_valid()) ++ordered_basis;
    }

public:
    OrderedBasis(const LieGroup& G) : ordered_basis{first_compatible_order(G)}, G{G} {
        advance_until_valid();
    }
    int hat(int i) const {
        return G.Dimension()-i-1;
    }
    OrderedBasis& operator++() {
        while (++ordered_basis && !is_valid()) {};
        return *this;
    } 
    exvector e() const {
        exvector e;
        std::transform(ordered_basis->begin(),ordered_basis->end(),back_inserter(e), [this] (int i) {return G.e()[i];} );
        return e;
    }
    operator bool() const {return ordered_basis;}
    vector<int> indices() const {return *ordered_basis;}
};


class Filtration {
    const LieGroup& G;
    OrderedBasis ordered_basis;    
    exvector w;    
    AlternativeLinearInequalities<WParameter> inequalities; 
    lst weights_;
/*
    AlternativeLinearInequalities<WParameter>(int i) const {
        ex sum=w_i+w_ihat - w_n;
        AlternativeLinearInequalities<WParameter>::nonzero(sum)
        LinearInequalities{}.positive(w_i+w_ihat - w_n) || LinearInequalities{}.negative(w_i+w_ihat - w_n)
    }

    multiplicities() const {
        for each i<n/2:

        auto alternatives = {w_i+w_ihat - w_n,w_n-w_i-w_ihat};
            w_i+w_ihat > w_n 
            w_iw_ihat == w_n, w_i<w_ihat , w_i=w_i-1
        
    }
*/
    AlternativeLinearInequalities<WParameter> multiplicity_at_least(int i, int multiplicity) const {
        AlternativeLinearInequalities<WParameter> result{};
        int delta=multiplicity-1;
        for (int k=max(i-delta,0); k<min(i+1,static_cast<int>(w.size())-delta);++k) {
            exvector x;
            for (int n=k;n<k+delta;++n) x.push_back(w[n]-w[n+1]);
            result = result || LinearInequalities<WParameter>{}.zero(x.begin(),x.end());
        }        
        return result;
    }

    AlternativeLinearInequalities<WParameter> inequalities_for_indices(int i, int ihat) const {
        assert(i!=ihat);
        auto result= AlternativeLinearInequalities<WParameter>::positive (w[i]+w[ihat]-w.back())
            || 
                AlternativeLinearInequalities<WParameter>::nonzero(w[i]-w[ihat]) &&
                    (multiplicity_at_least(i,2) || multiplicity_at_least(ihat,2))
            ||                
                AlternativeLinearInequalities<WParameter>::zero(w[i]-w[ihat]) && multiplicity_at_least(i,3);        
        return result;
    }

    AlternativeLinearInequalities<WParameter> inequalities_for_good_filtration() const {
        auto alt=inequalities_for_positive_filtration();    //filtration satisfying F1        
        exvector inequalities{}, strict_inequalities;
        for (int i=0;i<(w.size()+1)/2;++i) {           //F2
            inequalities.push_back(w[i]+w[ordered_basis.hat(i)]-w.back());        
            strict_inequalities.push_back(w[i]+w[ordered_basis.hat(i)]-w[w.size()-2]);
        }
        alt=alt &&  LinearInequalities<WParameter>{}.nonnegative(inequalities.begin(),inequalities.end()).positive(strict_inequalities.begin(),strict_inequalities.end());
        for (int i=0;i<(w.size())/2;++i)            //F3 and F4 for i!=ihat
            alt = alt &&inequalities_for_indices(i,ordered_basis.hat(i));        
        return alt;            
    }


    AlternativeLinearInequalities<WParameter> inequalities_for_positive_filtration() const {
        if (!ordered_basis) return {};
        exvector inequalities;        
        auto e=basis();       
        for (int i=1;i<e.size();++i)
            inequalities.push_back(w[i]-w[i-1]);        
        
        for (int i=0;i<G.Dimension();++i)
        for (int j=i+1;j<G.Dimension();++j)
        for (int k=0;k<G.Dimension();++k)
            if (!Hook(e[i]*e[j],G.d(e[k])).is_zero())
                inequalities.push_back(w[k]-w[i]-w[j]);
        return AlternativeLinearInequalities<WParameter>::positive(w[0]) &&        
            LinearInequalities<WParameter>{}.nonnegative(inequalities.begin(),inequalities.end());
    }
    void advance() {        
        if (++ordered_basis)
            inequalities=inequalities_for_good_filtration();
    }
    void advance_until_valid() {
        weights_=lst{};
        while (ordered_basis && (weights_=inequalities.find_solution()).nops()==0)
            advance();
    }
    static exvector create_parameters(int n) {
        exvector w;
        for (int i=1;i<=n;++i) w.push_back(WParameter(N.w(i)));
        return w;
    }
public:
    Filtration(const LieGroup& G) : G{G}, ordered_basis{G}, w{create_parameters(G.Dimension())}, inequalities{inequalities_for_good_filtration()} {              
        advance_until_valid();
    }

    const Filtration& operator++() {
        advance();
        advance_until_valid();
        return *this;
    }
    exvector basis() const {
        return ordered_basis.e();
    }    
    exvector weights() const {
        exvector result;
        transform(w.begin(),w.end(),back_inserter(result),[this] (ex w) {return w.subs(weights_);});
        //transform into integers by eliminating denominators
        auto take_lcm=[](ex x, ex y) {return lcm(x,y,false);};
        auto take_denominator=[](ex x) {return x.denom();};
        auto prod_of_denominators=transform_reduce(result.begin(),result.end(),ex{1},take_lcm,take_denominator);
        transform(result.begin(),result.end(),result.begin(),[prod_of_denominators] (ex x) {return x*prod_of_denominators;});
        return result;
    }
    operator bool() const {return ordered_basis;}
    vector<int> indices() const {return ordered_basis.indices();}
};

