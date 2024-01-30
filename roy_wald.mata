*! version 1.0.0 M.C. Park 24Jan2024

version 18
set matastrict on

mata:
void llku(transmorphic M, real rowvector b, real colvector lnf) {
real colvector p1, p2, p3, p4, p5, p6, p7, y, d

p1 = moptimize_util_xb(M, b, 1) // gamma
p2 = moptimize_util_xb(M, b, 2) // beta0
p3 = moptimize_util_xb(M, b, 3) // beta1
p4 = moptimize_util_xb(M, b, 4) // sig0
p5 = moptimize_util_xb(M, b, 5) // sig1
p6 = moptimize_util_xb(M, b, 6) // rho0
p7 = moptimize_util_xb(M, b, 7) // rho1

y  = moptimize_util_depvar(M, 1)
d  = moptimize_util_depvar(M, 2)

lnf = (1 :- d) :* (-ln(p4) :- .5 :* ((y :- p2) :/ p4):^2 :+ ///
lnnormal(-(p1 :+ (p6 :/ p4) :* (y :- p2)) :/ sqrt(1 :- p6:^2))) :+ ///
d :* (-ln(p5) :- .5 :* ((y :- p3) :/ p5):^2 :+ ///
lnnormal((p1 :+ (p7 :/ p5) :* (y :- p3)) :/ sqrt(1 :- p7:^2)))
}

class roy_wald {
	public:
	// initialize
	void new()
	
	// setup
	real scalar    N
	real scalar    kx
	real scalar    df
	real scalar    alp

	// estimates
	transmorphic   M
	void mleu()
	
	real colvector gam
	real colvector bet0
	real colvector bet1
	
	real scalar    sig0
	real scalar    sig1
	real scalar    rho0
	real scalar    rho1
	
	real matrix    Om
	
	// wald statistic
	real colvector coefs()
	real colvector h()
	real matrix    Dh()
	real scalar    chi2val()
	real scalar    pval()
	real scalar    reject()
	
	// report, returns
	void report()
	void returns()
}

void roy_wald::new() {
	N    = rows(st_data(., st_local("Y"), st_local("touse")))
	kx   = cols(st_data(., st_local("X"), st_local("touse"))) + 1
	df   = kx - st_numscalar("r(k_omitted)")
	alp  = .05
	M    = moptimize_init()
	mleu()
	
	gam  = moptimize_result_coefs(M)[1 .. kx]'
	bet0 = moptimize_result_coefs(M)[(kx + 1) .. (2 * kx)]'
	bet1 = moptimize_result_coefs(M)[(2 * kx + 1) .. (3 * kx)]'
	
	sig0 = moptimize_result_coefs(M)[(3 * kx + 1)]
	sig1 = moptimize_result_coefs(M)[(3 * kx + 2)]
	rho0 = moptimize_result_coefs(M)[(3 * kx + 3)]
	rho1 = moptimize_result_coefs(M)[(3 * kx + 4)]
	
	if (st_local("V") == "" | st_local("V") == "mle") {
		Om = moptimize_result_V(M)
	} else if (st_local("V") == "robust") {
		Om = moptimize_result_V_robust(M)
	}
}

void roy_wald::mleu() {
	moptimize_init_evaluator(M, &llku())
	
	moptimize_init_depvar(M, 1, st_data(., st_local("Y"), st_local("touse")))
	moptimize_init_depvar(M, 2, st_data(., st_local("D"), st_local("touse")))

	moptimize_init_eq_indepvars(M, 1, st_data(., st_local("X"), st_local("touse")))
	moptimize_init_eq_indepvars(M, 2, st_data(., st_local("X"), st_local("touse")))
	moptimize_init_eq_indepvars(M, 3, st_data(., st_local("X"), st_local("touse")))
	moptimize_init_eq_indepvars(M, 4, "")
	moptimize_init_eq_indepvars(M, 5, "")
	moptimize_init_eq_indepvars(M, 6, "")
	moptimize_init_eq_indepvars(M, 7, "")
	
	moptimize_init_tracelevel(M, "none")
	moptimize_init_valueid(M, "Log likelihood")
	
	moptimize_init_search(M, "on")
	moptimize_init_search_random(M, "off")
	moptimize(M)
}

real colvector roy_wald::coefs() {
	return((gam \ bet0 \ bet1 \ sig0 \ sig1 \ rho0 \ rho1))
}

real colvector roy_wald::h() {
	return(gam :* (sig1 * rho1 - sig0 * rho0) :- (bet1 :- bet0))
}

real matrix roy_wald::Dh() {
	return((I(kx) :* (sig1 * rho1 - sig0 * rho0), I(kx), -I(kx), ///
	-rho0 :* gam, rho1 :* gam, -sig0 :* gam, sig1 :* gam))
}

real scalar roy_wald::chi2val() {
	if (st_local("ADJ") == "" | st_local("ADJ") == "none") {
		return((h()' * invsym(Dh() * Om * Dh()') * h()))
	} else if (st_local("ADJ") == "df") {
		return(((N - length(coefs())) / N) * (h()' * invsym(Dh() * Om * Dh()') * h()))
	}
}

real scalar roy_wald::pval() {
	return(1 - chi2(df, chi2val()))
}

real scalar roy_wald::reject() {
	return(pval() < alp)
}

void roy_wald::report() {
	printf("\n")
	printf("{txt}%s \n", "{space 3}Test of the original Roy selection")
	printf("{txt}%s \n", "{space 3}Ho: selection rule is the original Roy")
	printf("\n")
	printf("{txt}%s {res}%g {txt}%s {res}%5.4f \n", "{space 3}chi2(", df, ")    = ", chi2val())
	printf("{txt}%s {res}%5.4f \n", "{space 3}Prob > chi2  = ", pval())
	if (reject() == 1) {
	printf("{txt}%s", "{space 3}Ho is rejected at the 5% level")
	} else if (reject() == 0) {
		printf("{txt}%s", "{space 3}Ho is not rejected at the 5% level")
	}
}

void roy_wald::returns() {
	st_numscalar("r(N)", N)
	st_numscalar("r(chi2)", chi2val())
	st_numscalar("r(df)", df)
	st_numscalar("r(p)", pval())
	st_numscalar("r(reject)", reject())
	
	if (st_local("V") == "") st_local("V", "mle")
	if (st_local("ADJ") == "") st_local("ADJ", "none")
}
end