/*20210513
*Specify the following in the calling dofile
*Variables of interest must be labeled (except outcome) and label values name must not match variable name or exceed 28 characters

*Specify outcome
local outcome trinfect

*Specify variables of interest
unab vars:  
disp "`vars'"

local first: word 1 of `vars'
disp "`first'"

local n_vars: word count `vars'
disp "`n_vars'"

*/

capture log close
log using "`results'/Analysis_`logpart'.log", replace

capture postutil clear

********************************************************************************

tempfile pvals
postfile test str35 varname str10 crude_chi2 str10 crude_pvalue using `pvals'

tempfile xtabs
postfile result str35 varname str35 value str35 overall str35 pos str35 crude_pr using `xtabs'

********************************************************************************
*Get overall prevalence
tab `outcome', matcell(Z)

*Local with total count
local tot_n = `r(N)'
di `tot_n'

*Local with value count if coded 0/1
local pos_n = Z[2,1]
di `pos_n'

*Local with percentage
local pos_p : di %3.2f `pos_n'/`tot_n'*100
di `pos_p'

if `den'==1 {
	local pos_den = " / " + "`tot_n'"
}
else {
	local pos_den = ""
}

local pos = "`pos_n'" + "`pos_den'" + " (" + "`pos_p'" + ")"

post result ("Overall") ("--") ("--") ("`pos'") ("--")

********************************************************************************
*Get prevalence by specified factors, crude prevalence ratios, and chi2 test results
local varct = 0
foreach var in `vars' {
	local varct = `varct' + 1

	use `analysis', clear

	local lblname : value label `var'
	di "`lblname' "

	tab `var' `outcome', matcell(X) matrow(Y) chi2 col

	local levelmax = `r(r)'

	*Local with total
	local tot_n = `r(N)'
	di `tot_n'
		
	local ivar = "i.`var'"
	
	poisson `outcome' `ivar', irr vce(cluster `vce')
	
	mat Z = r(table)

	testparm `ivar'
	
	local chi2 : di %3.2f `r(chi2)'
	
	if `r(p)'<0.01 {
		local pval = "<0.01"
	}
	else {
		local pval : di %3.2f `r(p)'
	}
	
	post test ("`var'") ("`chi2'") ("`pval'")
	
	forvalues i = 1/`levelmax' {

		*Store varname
		local var`i' = `"`var'"'
		di `" "`var`i''" "'

		*Scalar with current value (won't work with non-numeric)
		scalar val`i' = Y[`i',1]
		di val`i'

		local getval`i' = val`i'

		*Get value label for current value
		local vallbl`i' : label `lblname' `getval`i''
		di `" "`vallbl`i''" "'

		*Local with value total
		local val_n`i' = X[`i',1] + X[`i',2]
		di `val_n`i''
		
		*Local with value percentage of total
		local val_p`i' : di %3.2f `val_n`i''/`tot_n'*100
		di `val_p`i''
		
		*Local with value count - positive
		local pos_n`i' = X[`i',2]
		di `pos_n`i''

		*Local with pos denominator
		local pos_den`i' = `val_n`i''
		di `pos_den`i''
		
		*Local with percentage - positive
		local pos_p`i' : di %3.2f `pos_n`i''/`pos_den`i''*100
		di `pos_p`i''

		if `den'==1 {
			local overall_den`i' = " / " + "`tot_n'"
			local pos_den`i' = " / " + "`pos_den`i''"
		}
		else {
			local overall_den`i' = ""
			local pos_den`i' = ""
		}
		
		local overall`i' = "`val_n`i''" + "`overall_den`i''" + " (" + "`val_p`i''" + ")"
		
		local pos`i' = "`pos_n`i''" + "`pos_den`i''" + " (" + "`pos_p`i''" + ")"

		if `i'==1 {
			local pr`i' = "REF"
		}
		else {
			*Store coefficient and confidence interval
			local expb`i' : di %3.2f Z[1, `i']
			local ll`i' : di %3.2f Z[5, `i']
			local ul`i' : di %3.2f Z[6, `i']

			local pr`i' = "`expb`i''" + " (" + "`ll`i''" + ", " + "`ul`i''" + ")"
		}
		
		post result ("`var`i''") ("`vallbl`i''") ("`overall`i''") ("`pos`i''") ("`pr`i''")

	}
}
postclose result
postclose test

use `xtabs', clear
gen n = _n
merge m:1 varname using `pvals', gen(merge)
replace crude_pr = "--" if merge==1
replace crude_chi2 = "--" if merge==1
replace crude_pvalue = "--" if merge==1
drop merge

save `xtabs', replace

********************************************************************************
*Get final adjusted model based on all possible subsets

use `analysis', clear

*Make a macro called tuples with all of your control variables
tuples `vars', asis display

tempfile allsubsets
tempname memhold
postfile models str100 model obs aic bic using `allsubsets'

forvalues i=1/`ntuples'{ 
	quietly poisson nonresponse `tuple`i'', irr vce(cluster `vce') 
	quietly estat ic
	mat s=r(S)
	local obs=s[1,1]
	local aic=s[1,5]
	local bic=s[1,6]
	post models ("`tuple`i''") (`obs') (`aic') (`bic') 
}
postclose models

use `allsubsets', clear

*Identify final model based on minimum BIC
sort bic
local selvars = model[1]
disp "`selvars'"

********************************************************************************
*Get final model results
use `analysis', clear

local finalvars 
foreach var in `selvars' {
	local ivar = "i.`var'"
	local finalvars = "`finalvars'" + "`ivar' "
}
	
poisson `outcome' `finalvars', irr vce(cluster `vce') 
	
mat Z = r(table)
mat li Z

tempfile finalpvals
postfile test2 str35 varname str10 final_chi2 str10 final_pvalue using `finalpvals'

foreach var in `finalvars' {
	testparm `var'
	
	local varname = subinstr("`var'", "i.", "", 1)
	
	local chi2 : di %3.2f `r(chi2)'
	
	if `r(p)'<0.01 {
		local pval = "<0.01"
	}
	else {
		local pval : di %3.2f `r(p)'
	}
	
	post test2 ("`varname'") ("`chi2'") ("`pval'")
}

tempfile finalresults
postfile result2 str35 varname str35 value str35 final_pr using `finalresults'

local colsZ : colsof Z
local colsZnocons = `colsZ'-1
local colnamesZ : colnames Z
forvalues j = 1/`colsZnocons' {
	
	local colname : word `j' of `colnamesZ'
	di "`colname'"
	local colval = substr("`colname'", 1, 1)
	di "`colval'"
	local pos = strpos("`colname'", ".") +1
	local varname = substr("`colname'", `pos', .)
	di "`varname'"
	*Get label name for variable
	local lblname : value label `varname'
	di "`lblname'"
	*Get value label for current value
	local vallbl : label `lblname' `colval'
	di `" "`vallbl'" "'
	
	if `colval'==1 {
		local pr = "REF"
	}
	else {
		*Store coefficient and confidence interval
		local expb : di %3.2f Z[1, `j']
		local ll : di %3.2f Z[5, `j']
		local ul : di %3.2f Z[6, `j']

		local pr = "`expb'" + " (" + "`ll'" + ", " + "`ul'" + ")"
	}

	post result2 ("`varname'") ("`vallbl'") ("`pr'")
}
postclose test2
postclose result2

use `finalresults', clear
merge m:1 varname using `finalpvals', nogen

save `finalresults', replace

use `xtabs', clear
merge 1:1 varname value using `finalresults', gen(finalmerge)

log close

********************************************************************************
*Format final table output

replace final_pr = "--" if finalmerge==1
replace final_chi2 = "--" if finalmerge==1
replace final_pvalue = "--" if finalmerge==1

drop finalmerge

sort n

bysort varname (n): gen first = _n == 1

expand 2 if first

sort n

gen n2 = _n

bysort var (n2): gen n3 = _n

foreach var of varlist value overall pos crude_pr final_pr {
	replace `var' = "" if n3==1

}
foreach var of varlist varname crude_chi2 crude_pvalue final_chi2 final_pvalue {
	replace `var' = "" if n3!=1
}

sort n2

drop first n* 

compress
