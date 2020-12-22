/*===========================================================================================*/
/*HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH*/
/*===========================================================================================*/
/*						Code for statistical frontier analysis	 	 	 	 	 	 	 	 */
/*===========================================================================================*/
/*HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH*/
/*===========================================================================================*/

// Data
use coeff_by_model, clear

// Dependent variable: Additional fuel consumption due to higher temperatures (L per 100 km)
// Co-variates: Displacement (L), Maximum horsepower (1,000 kW), Weight (ton), Wheelbase (m), Type-approval fuel efficiency (L/100 km), Market price (10,000 yuan)
keep coeff_hot displacement weight ref_price ref_fe wheelbase maxhorsepower _freq trip_distance

capture drop uci_exp uci_exp_LB95 uci_exp_UB95 uci_exp2 uci_exp2_LB95 uci_exp2_UB95 uci_exp3 uci_exp3_LB95 uci_exp3_UB95
gen a_weight = _freq * trip_distance


/*==================================================================================*/
/*				Basic model - Exponential inefficiency term							*/
/*==================================================================================*/

// All co-variates
sfcross coeff_hot displacement ref_price ref_fe maxhorsepower wheelbase weight, distribution(exponential) cost nolog robust
eststo m1
predict uci_exp1, u ci
reg uci_exp1 [aweight = a_weight]

/*==================================================================================*/
/*				Basic model - Half normal inefficiency term							*/
/*==================================================================================*/

// All co-variates
sfcross coeff_hot displacement ref_price ref_fe maxhorsepower wheelbase weight, distribution(hnormal) cost nolog robust
eststo m2
predict uci_exp2, u ci
reg uci_exp2 [aweight = a_weight]


// Two key co-variates
sfcross coeff_hot ref_price maxhorsepower, distribution(hnormal) cost nolog robust
eststo m3
predict uci_exp3, u ci
reg uci_exp3 [aweight = a_weight]

/*==================================================================================*/
/*			Model with heteroskedasticity - Half normal inefficiency term			*/
/*==================================================================================*/

sfcross coeff_hot displacement ref_price ref_fe maxhorsepower wheelbase weight, distribution(hnormal) vsigma(ref_price) cost nolog robust
eststo m4
matrix vec_sfcross = e(b)
matrix vec_usigma0 = vec_sfcross[1, "Usigma:"]
predict uci_exp4, u ci
reg uci_exp4 [aweight = a_weight]

/*==================================================================================*/
/*			Rebound-effects model - Half normal inefficiency term					*/
/*==================================================================================*/

sfcross coeff_hot displacement ref_price ref_fe maxhorsepower wheelbase weight, distribution(hnormal) usigma(ref_price) vsigma(ref_price) cost nolog robust
eststo m5
predict pred_i, xb
gen epsilon_i = coeff_hot - pred_i

// Save estimates
matrix vec_sfcross_h = e(b)
matrix vec_frontier = vec_sfcross_h[1, "Frontier:"]
matrix vec_usigma_cons = vec_sfcross_h[1, "Usigma:_cons"]
scalar usigma_cons = vec_usigma_cons[1,1]
matrix vec_usigma_ref_price = vec_sfcross_h[1, "Usigma:ref_price"]
scalar usigma_ref_price = vec_usigma_ref_price[1,1]
matrix vec_vsigma_cons = vec_sfcross_h[1, "Vsigma:_cons"]
scalar vsigma_cons = vec_vsigma_cons[1,1]
matrix vec_vsigma_ref_price = vec_sfcross_h[1, "Vsigma:ref_price"]
scalar vsigma_ref_price = vec_vsigma_ref_price[1,1]

egen usigma0 = max(usigma_cons  + usigma_ref_price * ref_price)
gen scaling = sqrt(exp((usigma_cons  + usigma_ref_price * ref_price) - usigma0))	// Rebound effect
gen exp_sigma_u0 = sqrt(exp(usigma0))
gen sigma_u = usigma_cons  + usigma_ref_price * ref_price
gen exp_sigma_u = sqrt(exp(usigma_cons + usigma_ref_price * ref_price))
gen sigma_v = vsigma_cons  + vsigma_ref_price * ref_price
gen exp_sigma_v = sqrt(exp(vsigma_cons + vsigma_ref_price * ref_price))
gen exp_sigma_i = sqrt(exp_sigma_v * exp_sigma_v + exp_sigma_u * exp_sigma_u)

gen u_i = epsilon_i * (exp_sigma_u * exp_sigma_u) / (exp_sigma_i * exp_sigma_i)
gen sigma_i = exp_sigma_u * exp_sigma_v / exp_sigma_i

gen E_ui_tilt = u_i + sigma_i * (normalden(-u_i/sigma_i) / (1 - normal(-u_i/sigma_i)))	// Overall inefficiency
reg E_ui_tilt [aweight = a_weight]

esttab using sfa_results.csv, b(3) se(2) replace	// Table SXX
eststo clear

// Show correlation of estimates across models
corr uci_exp2 uci_exp3 uci_exp4 E_ui_tilt



