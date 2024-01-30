*! version 1.0.0 M.C. Park 24Jan2024

program define roy_wald, rclass byable(recall)
	version 18
	
	syntax varlist(min = 2 numeric fv) [if] [in], SELect(varname) [vce(string) adj(string)]
	marksample touse
	markout `touse' `select'
	
	gettoken Y X : varlist
	gettoken D   : select
	gettoken V   : vce
	gettoken ADJ : adj
	
	_rmcoll `X', expand // omitted variables check

	mata: instance = roy_wald()
	mata: instance.report()
	mata: instance.returns()

	return scalar N      = `r(N)'
	return scalar chi2   = `r(chi2)'
	return scalar df     = `r(df)'
	return scalar p      = `r(p)'
	return scalar reject = `r(reject)'
	
	return local  vce    = "`V'"
	return local  adj    = "`ADJ'"
end

include roy_wald.mata, adopath