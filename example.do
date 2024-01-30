/* roy_wald example: marriage choice */
webuse nlswork, clear
generate agesq = age^2

// data description
describe
summarize

// basic usage (pooled sample)
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours, sel(msp)
return list

// robust vce
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours, sel(msp) vce(robust)

// adjustment for degrees of freedom
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours, sel(msp) adj(df)

// subsample (in)
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours in 1/3000, sel(msp) vce(robust)

// subsample (if)
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours if year == 68, sel(msp) vce(robust)
roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours if year == 88, sel(msp) vce(robust)

// byable
keep if year == 68 | year == 88
bysort year: roy_wald ln_wage age agesq ib3.race collgrad ttl_exp hours, sel(msp) vce(robust)

// help file
help roy_wald