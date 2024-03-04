
//iterates through linear extensions of a poset in lexicographic order wrt some other total ordering <
//the partial ordering is passed in the constructor as a subset of the Cartesian product
//the total ordering to use for lexicographic ordering is T::operator<(T)
template<typename T>
struct LinearExtension {
    set<T> to_add; //use of set ensures that iteration through to_add is done relative to the total ordering. public functions leave this empty
    set<pair<T,T>> partial_order_relation;  
    vector<T> order;

    bool less_than_in_partial_order(T first, T second) const {
        return partial_order_relation.count(make_pair(first,second));
    }

    bool has_incoming_edges_from_to_add(T node) const {
        for (auto n: to_add)
            if (less_than_in_partial_order(n,node)) return true;
        return false;
    }
    T first_element_with_no_incoming_edges_from_to_add() const {
        for (auto n: to_add)
            if (!has_incoming_edges_from_to_add(n)) return n;
        throw logic_error("error in LinearExtension: not_a_poset");
    }
    void complete() {
        while (!to_add.empty()) {
            auto node=first_element_with_no_incoming_edges_from_to_add();
            to_add.erase(node);
            order.push_back(node);
        }
    }
    optional<T> first_element_with_no_incoming_edges_from_to_add_and_greater_than(T node) const {
        for (auto n: to_add)
            if (n>node && !has_incoming_edges_from_to_add(n)) return n;
        return nullopt;
    }    
    LinearExtension(const vector<T>& poset, const set<pair<T,T>>& partial_order_relation) : to_add(poset.begin(),poset.end()), partial_order_relation{partial_order_relation} {
        complete();
    }
    LinearExtension()=default;
public:
    static LinearExtension begin(const vector<T>& poset, const set<pair<T,T>>& partial_order_relation) {
        return LinearExtension{poset,partial_order_relation};
    }
    static LinearExtension end(const vector<T>& poset, const set<pair<T,T>>& partial_order_relation) {
        return LinearExtension{};
    }


//take the next linear extension in lexicographic order. 
    LinearExtension& operator++() {
        //starting from the end, remove elements from order, and try to replace them with higher values
        optional<T> next;
        while (!next && !order.empty()) {
            auto last=order.back();
            order.pop_back();
            to_add.insert(last);
            next=first_element_with_no_incoming_edges_from_to_add_and_greater_than(last);
        }
        if (next) {
            to_add.erase(next.value());
            order.push_back(next.value());
            complete();        
        }        
        return *this;
    }
    bool operator!=(const LinearExtension& other) const  {
        return order!=other.order;
    }
    const vector<T>& operator*() const {return order;}
    const vector<T>* operator->() const {return &order;}
    operator bool() const {return !order.empty();}
};



template<typename T>
class LinearExtensions {
    LinearExtension<T> begin_,end_;
public:
    LinearExtensions(const vector<T>& poset, const set<pair<T,T>>& partial_order_relation) 
        : begin_{LinearExtension<T>::begin(poset,partial_order_relation)}, end_{LinearExtension<T>::end(poset,partial_order_relation)} {}
    auto begin() const {return begin_;}
    auto end() const {return end_;}
};