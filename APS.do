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

*capture log close
*log using "`results'/Analysis_`logpart'.log", replace

capture postutil clear

********************************************************************************

tempfile xtabs
postfile result str35 varname str35 valuelabel str35 value str35 overall str35 pos str35 expb str10 pvalue str7 tablepart str2 n2 using `xtabs'

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

post result ("Overall") ("--") ("--") ("--") ("`pos'") ("--") ("--") ("overall") ("")

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
	
	post result ("`var'") ("") ("") ("") ("") ("") ("`pval'") ("stat") ("")
	
	forvalues i = 1/`levelmax' {

		*Store varname
		local var = `"`var'"'
		di `" "`var'" "'

		*Current level value
		local value = Y[`i',1]
		di `" "`value'" "'

		*Get value label for current value
		local vallbl : label `lblname' `value'
		di `" "`vallbl'" "'

		*Local with row marginal total
		local val_n = X[`i',1] + X[`i',2]
		di `val_n'
		
		*Local with row percentage of overall total
		local val_p : di %3.2f `val_n'/`tot_n'*100
		di `val_p'
		
		*Local with count of positives within level
		local pos_n = X[`i',2]
		di `pos_n'

		*Local with pos denominator, uses row marginal total
		local pos_den = `val_n'
		di `pos_den'
		
		*Local with percentage - positive
		local pos_p : di %3.2f `pos_n'/`pos_den'*100
		di `pos_p'

		if `den'==1 {
			local overall_den = " / " + "`tot_n'"
			local pos_den = " / " + "`pos_den'"
		}
		else {
			local overall_den = ""
			local pos_den = ""
		}
		
		local overall = "`val_n'" + "`overall_den'" + " (" + "`val_p'" + ")"
		
		local pos = "`pos_n'" + "`pos_den'" + " (" + "`pos_p'" + ")"

		if `i'==1 {
			local pr = "REF"
		}
		else {
			*Store coefficient and confidence interval
			local expb : di %3.2f Z[1, `i']
			local ll : di %3.2f Z[5, `i']
			local ul : di %3.2f Z[6, `i']

			local pr = "`expb'" + " (" + "`ll'" + ", " + "`ul'" + ")"
		}
		
		post result ("`var'") ("`vallbl'") ("`value'") ("`overall'") ("`pos'") ("`pr'") ("") ("model") ("`value'")
	}
} 
postclose result

use `xtabs', clear
/********************************************************************************
*Get final adjusted model based on all possible subsets

use `analysis', clear

*Make a macro called tuples with all of your control variables
tuples `vars', asis //display

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

********************************************************************************/
*Get final model results
use `analysis', clear

*TEMPORARY STEP
local selvars = "ses site_id mobile_phone"

local finalvars 
foreach var in `selvars' {
	local finalvars = "`finalvars'" + " i.`var'"
}

poisson `outcome' `finalvars', irr vce(cluster `vce') 
	
*Output matrix with regression results
mat Z = r(table)'
mat li Z

tokenize "`:rownames Z'" //tokenize splits string into pieces referred to by `1',`2', etc
local rows : rowsof Z //get count of rows

preserve

clear
svmat Z, names(col) //convert matrix to dataset
gen rowname = ""

forvalues t=1/`rows' {
	replace rowname = "``t''" if _n==`t' //populate variable with rownames extracted from matrix
}

drop if rowname=="_cons"

split rowname, parse(".")
rename rowname2 varname
	
tostring b, replace format(%3.2f) force
tostring ll, replace format(%3.2f) force
tostring ul, replace format(%3.2f) force

gen final_expb = "REF"
replace final_expb = b + " (" + ll + ", " + ul + ")" if !strpos(rowname1,"b")

replace rowname1 = subinstr(rowname1,"b","",.)
replace rowname1 = subinstr(rowname1,"o","",.)
rename rowname1 n2
		
keep varname final_expb n2

tempfile final
save `final'

restore

tempfile FINAL_p
postfile FINAL_p str35 varname str10 final_p str7 tablepart using `FINAL_p'

foreach var in `selvars' {
	
	di "Test `var'"
	testparm i.`var'
			
	if `r(p)'<0.01 {
		local pval = "<0.01"
	}
	else {
		local pval : di %3.2f `r(p)'
	}
	
	post FINAL_p ("`var'")  ("`pval'") ("stat")
}
postclose FINAL_p

********************************************************************************
if `log'==1 {
log off
}
*Format final table output
use `xtabs', clear

gen n = _n

preserve

keep if tablepart!="model"

merge 1:1 varname tablepart using `FINAL_p', gen(finalpmerge)

tempfile other
save `other'

restore

keep if tablepart=="model"

merge 1:1 varname n2 using `final', gen(finalmerge)

append using `other'

sort n

drop finalpmerge finalmerge

bysort varname (n2): replace varname = "" if _n!=1

sort n

drop tablepart n n2

order varname valuelabel value overall pos expb pvalue final_expb final_p

compress

save `xtabs', replace

if `log'==1 {
log on
}
