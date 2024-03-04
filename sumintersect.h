Subspace<DifferentialForm> intersected_with(const VSpace<DifferentialForm>& W, const VSpace<DifferentialForm>& V) {
	auto V_in_sum=Subspace<DifferentialForm>(V.e_begin(),V.e_end(),W.e_begin(),W.e_end());
	auto w=W.GenericElement();
	lst eqns;
	GetCoefficients<DifferentialForm>(eqns,V_in_sum.ProjectOnComplement(w));
	return W.SubspaceFromEquations(eqns.begin(),eqns.end());
}

VectorSpace<DifferentialForm> sum(const VSpace<DifferentialForm>& V, const VSpace<DifferentialForm>& W) {
	exvector all_generators;
	all_generators.insert(all_generators.end(),V.e_begin(),V.e_end());
	all_generators.insert(all_generators.end(),W.e_begin(),W.e_end());
	return {all_generators};
}
