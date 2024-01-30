{smcl}
{* *! version 1.0.0 30Jan2021}
{viewerjumpto "options" "roy_wald##options"}{...}
{title:Title}

{phang}
{bf:roy_wald} {hline 1} A test for the Roy selection mechanism

{title:Syntax}

{phang}
{cmd:roy_wald} {depvar} {indepvars} {ifin} , {opth sel:ect(varname)} [{cmd:}{help roy_wald##options:options}]
{p_end}

{synoptset 16 tabbed}{...}
{synopt: {it:depvar}} dependent variable in outcome equations{p_end}
{synopt: {it:indepvars}} independent variables in outcome equations{p_end}
{synopt: {it:varname}} dependent variable in the selection equation{p_end}

{title:Description}

{pstd}
{bf:roy_wald} tests the Roy selection mechanism. In particular, the null hypothesis is the original Roy selection rule, and the alternative is the extended or generalized Roy selection. Non-rejection of the null hypothesis indicates that the outcome variable is the only determinant of agents' choice.


{marker options}{...}
{title:Options}

{dlgtab:Covariance Matrix}

{phang}
{opt vce}({it:vcetype}) specifies the type of the covariance matrix used in the Wald statistic, which includes the standard ({opt mle}) and the sandwich ({opt robust)} covariance matrices. The default type is {opt mle}.
{p_end}

{dlgtab:Adjustment for Degrees of Freedom}

{phang}
{opt adj}({it:string}) specifies the small-sample adjustment. The default is to not use the adjustment ({opt none}), and the adjustment for degrees of freedom ({opt df}) is available.
{p_end}


{title:Examples}

	1) Basic usage

	{cmd:roy_wald Y X1 X2 ... Xk, select(D)}

	2) Robust covariance matrix

	{cmd:roy_wald Y X1 X2 ... Xk, select(D) vce(robust)}

	3) Adjustment for degrees of freedom

	{cmd:roy_wald Y X1 X2 ... Xk, select(D) adj(df)}

	4) By group

	{cmd:bysort varname: roy_wald Y X1 X2 ... Xk, select(D)}


{title:Author}

{phang}
Minchul Park (Korea University Economics), minchul1352@korea.ac.kr
{p_end}