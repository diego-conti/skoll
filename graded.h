#include <set>

class empty_sequence_t {} empty_sequence;

struct comes_before_in_list {
    exvector list;
    bool operator()(ex x, ex y) const {
        return find(list.begin(),list.end(),x)<find(list.begin(),list.end(),y);
    }  
};


class PartialWeightSequence {
    using Compare = ex_is_less;
    exvector w;
    //comes_before_in_list comparator;
    Compare comparator;
    multiset<ex,Compare> unassigned_weights;
    int multeplicity(ex wi) const {
        return count(w.begin(),w.end(),wi)+count(unassigned_weights.begin(),unassigned_weights.end(),wi);
    }

    //return true if w_j is a weight and adding w_k would violate (G1), since w_i+w_j=w_k violates (G1) and i is not less than k
    bool violates_G1(ex w_i, ex w_j, ex w_k) const {
        if (w_i==w_j)
            return multeplicity(w_i)>1;
        else return multeplicity(w_j);
    }

    bool adding_preserves_invariance(ex w_k) const {
        for (auto w_i : unassigned_weights)
            if (violates_G1(w_i,w_k-w_i,w_k)) return false;
        return true;
    }

    bool valid() const {return !w.empty();}

    exvector with_added(ex w_k) const {
        PartialWeightSequence result(*this);
        result.w.push_back(w_k);
        result.unassigned_weights.extract(w_k);
        return result.first_complete();        
    }
    ex remove_last() {        
        auto w_k=w.back();
        w.pop_back();
        unassigned_weights.insert(w_k);
        return w_k;
    }

public:
    PartialWeightSequence(const exvector& weights) : 
        w{weights} {
    }
    PartialWeightSequence(const exvector& weights,empty_sequence_t) :// comparator{weights},
        unassigned_weights{weights.begin(),weights.end(),comparator}    {
            assert(!weights.empty());
    }
    exvector first_complete() const {
        if (unassigned_weights.empty()) return w;
        for (auto w_k: unassigned_weights) {
            //cout<<horizontal(w)<<" "<<w_k<<endl;
            if (adding_preserves_invariance(w_k)) {
                auto with_k=with_added(w_k);                
              //  cout<<horizontal(with_k)<<endl;
                if (!with_k.empty()) return with_k;                
            }
        }
        return {};
    }
    exvector next_complete(ex lbound_for_wk) const {
        for (auto w_k: unassigned_weights) {
            //cout<<"completing: "<<horizontal(w)<<" "<<w_k<<endl;
            if (comparator(lbound_for_wk,w_k) && adding_preserves_invariance(w_k)) {            
                auto with_k=with_added(w_k);   
              //  cout<<horizontal(with_k)<<endl;             
                if (!with_k.empty()) return with_k;
                else lbound_for_wk=w_k; //avoid trying duplicates (taking advantage of the fact that the multiset is ordered by comparator)
            }
        }
        return {};
    }
    exvector next() {
        exvector next;
        while (!w.empty() && next.empty()) {
            auto last=remove_last();
            //cout<<"removed "<<horizontal(w)<<" "<<last<<endl;
            next=next_complete(last);                        
        }
        return next;
    }
};

class Iterator {
    void advance_until_valid() {
        do {
            PartialWeightSequence p{w};
            w=p.next();
        }
        while (!w.empty() && !is_valid());
    }
    bool wi_plus_wihat_equals_wj_violates(int i, int ihat, int j) const { //if w[i]+w[ihat]=w[j]
        if (j<w.size()-1) return true; //violates (G2)
        if (w[i]!=w[ihat] && multiplicity(w[i])==1 && multiplicity(w[ihat])==1) return true; //violates (G3)
        if (w[i]==w[ihat] && i!=ihat && multiplicity(w[i])<=2) return true; //violates (G4)
        return false;
    }

    bool is_valid() const {
        for (int i=0;i<(w.size()+1)/2;++i) {
            auto ihat=w.size()-i-1;
            auto wi_plus_wihat=w[i]+w[ihat];
            auto j=find(w.begin(),w.end(),wi_plus_wihat)-w.begin();
            if (j<w.size() && wi_plus_wihat_equals_wj_violates(i,ihat,j)) return false;            
        }
        for (int i=0;i<w.size();++i) {
            auto ihat=w.size()-i-1;
            for (int j=0;j<(w.size()+1)/2;++j) {
                auto jhat=w.size()-j-1;
                if (multiplicity(w[i]+w[j]) && multiplicity(w[ihat]+w[jhat]) && (j!=ihat || i!=jhat)) return false; //violates (G5)
            }
        }
        return true;
    }
    exvector w;
    int multiplicity(ex weight) const {
        return count(w.begin(),w.end(),weight);
    }
public:
    static Iterator begin(const exvector& weights) {
        Iterator result;
        PartialWeightSequence p{weights,empty_sequence};
        result.w=p.first_complete();
        if (result.w.size() && !result.is_valid()) result.advance_until_valid();
        return result;
    }
    static Iterator end(const exvector& weights) {
        return Iterator{};
    }
    Iterator& operator++() {
        advance_until_valid();
        return *this;
    }    
    bool operator!=(const Iterator& o) const {
        return w!=o.w;
    }
    Iterator operator++(int) {
        Iterator result(*this);
        ++*this;
        return result;
    }
    exvector operator*() const {return w;}
};

class WeightSequencesRespectingOrder {
    exvector weights;
public:
    WeightSequencesRespectingOrder(const exvector& weights) : weights{weights} {     
    }
    Iterator begin() const {
        return (weights.empty())? Iterator::end(weights) : Iterator::begin(weights);
    }
    Iterator end() const {
        return Iterator::end(weights);
    }
    vector<int> as_indices(const exvector& weights) const {
        auto original_weights=this->weights;
        auto find_in_weights=[&original_weights] (ex w){
            auto i=find(original_weights.begin(),original_weights.end(),w);
            *i=lst();   //empty placeholder which equals no weight
            return i-original_weights.begin();
        };
        vector<int> result;
        transform(weights.begin(),weights.end(),back_inserter(result),find_in_weights);
        return result;
    }
};
