struct Couple {
	int first;
	int second;
	operator bool() const {return first!=second;}
	Couple& advance(int n) {
		if (second<n-1) ++second;
		else if (first<n-2) {
			++first;
			second=first+1;
		}
		else first=second=0;
		return *this;
	}
	bool overlaps(Couple x) {
		return first==x.first || second==x.first||first==x.second || second==x.second;
	}
	bool operator!=(const Couple& other) const {
		return first!=other.first || second!=other.second;
	}
};

class CoupleList {
	int n;	
	list<Couple> couples;
	bool compatible_from(list<Couple>::iterator i) {
		for (;i!=couples.end();++i) {
			auto j=i;
			while (++j!=couples.end())
				if (i->overlaps(*j)) return false;
		}
		return true;
	}
	
	bool advance_till_compatible_to_the_right(list<Couple>::iterator i) {
		do {
			i->advance(n);
		}
		while (*i && !compatible_from(i));
		return *i;	
	}
	bool advance_till_compatible(list<Couple>::iterator i) {
		if (!advance_till_compatible_to_the_right(i)) return false;
		for (auto j=i;j!=couples.begin();) {
			*--j=*i;
			if (!advance_till_compatible_to_the_right(j)) return false;
		}
		return true;
	}
public:
	CoupleList() : n{0} {};
	CoupleList(int n, int k) : n{n} {
		assert(n>=2*k);
		while (k--) {
			couples.push_front({0,1});
			if (!compatible_from(couples.begin())) advance_till_compatible(couples.begin());
		}		
	}
	CoupleList(int n,std::initializer_list<Couple> list) : n{n},couples{list} {}
	CoupleList& operator++() {
		auto iter=couples.begin();	
		assert(iter!=couples.end());
		while (!advance_till_compatible(iter)) {
			if (++iter==couples.end()) {
				couples.clear();
				break;
			}
		}
		return *this;
	}
	auto begin() const {return couples.begin();}
	auto end() const {return couples.end();}
	bool empty() const {return couples.empty();}	
	bool operator!=(const CoupleList& other) const {		
		return couples.size()!=other.couples.size() || 
			std::any_of(couples.begin(),couples.end(), [other] (auto couple) {
				return find(other.begin(),other.end(),couple)==other.end();
			});
	}
};

class OrderTwoAutomorphism {
	int n;
	CoupleList couples;	
	static int apply(int k, Couple trasposition) {
		if (k==trasposition.first) k= trasposition.second;
		else if (k==trasposition.second) k= trasposition.first;
		return k;
	}
public:
	OrderTwoAutomorphism() : n{0} {};
	OrderTwoAutomorphism(int n, int k) : n{n}, couples{n,k} {
	}
	OrderTwoAutomorphism(int n,std::initializer_list<Couple> list) : n{n},couples{n,list} {}
	OrderTwoAutomorphism& operator++() {
		++couples;
		return *this;
	}
	operator bool() const {return !couples.empty();}
	bool operator!=(const OrderTwoAutomorphism& other) const {		
		return n!=other.n || couples!=other.couples;
	}	
	int apply(int node) const {
		for (auto& couple : couples) node=apply(node,couple);
		return node;	
	}
	string to_string() const {
		string result;
		for (auto c: couples) result+="("+std::to_string(c.first+1)+" "+std::to_string(c.second+1)+")";
		if (couples.empty()) result="()";
		return result;	
	}
	matrix to_matrix() const {
		matrix sigma(n,n);
		for (int i=0;i<n;++i) sigma(i,i)=1;
		for (auto couple: couples) {
			sigma(couple.first,couple.second)= sigma(couple.second,couple.first)=1;
			sigma(couple.first,couple.first)= sigma(couple.second,couple.second)=0;
		}
		return sigma;
	}
	matrix sigma_diagonal_metric(const exvector& coefficients) const {
		matrix sigma(n,n);
		for (int i=0;i<n;++i) sigma(i,i)=coefficients[i];
		for (auto couple: couples) {
			sigma(couple.first,couple.second)= sigma(couple.second,couple.first)=coefficients[couple.first];
			sigma(couple.first,couple.first)= sigma(couple.second,couple.second)=0;
		}
		return sigma;
	}	
};

inline ostream& operator<<(ostream& os, const OrderTwoAutomorphism& sigma) {
	return os<<sigma.to_string();
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(MetricParameter,symbol)

class SigmaDiagonalMetricIterator {
	int n,k;
	OrderTwoAutomorphism sigma;
	exvector coefficients;
public:
	SigmaDiagonalMetricIterator& operator++() {
		if (!sigma && k==0) {
			if (++k<=n/2) sigma=OrderTwoAutomorphism{n,k};
		}
		else {
			assert(sigma);
			++sigma;		
			if (!sigma && ++k<=n/2) sigma=OrderTwoAutomorphism{n,k};
		}
		return *this;
	}
	matrix operator*() const {return sigma.sigma_diagonal_metric(coefficients);}	

	static SigmaDiagonalMetricIterator begin(int n) {
		SigmaDiagonalMetricIterator result;	
		result.n=n; result.k=0;	result.sigma=OrderTwoAutomorphism{n,0};
		for (int i=1;i<=n;++i) result.coefficients.push_back(MetricParameter(N.g(i)));
		return result;		
	}
	static SigmaDiagonalMetricIterator end(int n) {
		SigmaDiagonalMetricIterator result;	
		result.n=n; result.k=n/2+1;	 result.sigma=OrderTwoAutomorphism{n,0};
		return result;
	}
	bool operator!=(const SigmaDiagonalMetricIterator& other) const {		
		return n!=other.n || k!=other.k || sigma!=other.sigma;
	}	
	string to_string() const {
		stringstream s;
		s<<n<<": "<<sigma.to_string();
		return s.str();
	}
};


class SigmaDiagonalMetrics {
	int n;
public:
	SigmaDiagonalMetrics(int n) : n{n} {}
	auto begin() const {
		return SigmaDiagonalMetricIterator::begin(n);
	}
	auto end() const {
		return SigmaDiagonalMetricIterator::end(n);
	}
};