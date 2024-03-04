WEDGE_DECLARE_NAMED_ALGEBRAIC(Parameter,symbol)

class GlobalParameters {
    exvector parameters;
public:
    GlobalParameters(const Name& name, int n=100) {
        for (int i=0;i<n;++i)
            parameters.push_back(Parameter(name(i)));            
    }
    ex operator[](int i) const {return parameters[i];}
};
matrix antidiagonal_metric(const vector<int>& indices)  {
    static GlobalParameters parameters(N.g);
    matrix g(indices.size(),indices.size());    
    auto i=indices.begin();
    auto j=indices.rbegin();
    for (; i!=indices.end();++i,++j)
        g(*i,*j)=g(*j,*i)=parameters[*i+1];
    return g;
}