/*===========================================================================================*/
/*HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH*/
/*===========================================================================================*/
/*				Code for main regressions and robustness checks	 	 	 	 	 	 	 	 */
/*===========================================================================================*/
/*HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH*/
/*===========================================================================================*/

/*==================================================================================*/
/*		        	Regressions: impact of temp	on fuel economy     				*/
/*==================================================================================*/
// Data
use data, clear
xtset uuid edate_hour

// Examine the distribution of temperature
sum avg_avg_temp, detail
histogram avg_avg_temp, bin(20) fraction xlabel(#10)

// The first set of temperature bins: 10 bins   //the 5th bin (15-20˚C), coded as 4, is the baseline group
egen avg_avg_temp_id = cut(avg_avg_temp), at(-30,0,5,10,15,20,25,27,29,31,40) icodes
tab avg_avg_temp_id

xtreg est_fe i(0/3 5/9)bn.avg_avg_temp_id, fe robust cluster(brand_group)
eststo m1
xtreg est_fe i(0/3 5/9)bn.avg_avg_temp_id i.year_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo m2
xtreg est_fe i(0/3 5/9)bn.avg_avg_temp_id i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo m3

estadd local Year_fixed_effects "NO":m1
estadd local Year_fixed_effects "YES":m2 m3
estadd local Month_fixed_effects "NO":m1 m2
estadd local Month_fixed_effects "YES":m3
estadd local Vehicle_fixed_effects "YES": m1 m2 m3

esttab m1 m2 m3 using ac_table1.csv, label b(3)  se(3) ar2 drop(*.month_mid *.year_mid) nocon ///
		compress star(* 0.1 ** 0.05 *** 0.01) ///
		scalars(Year_fixed_effects Month_fixed_effects Vehicle_fixed_effects)

// Get the weighting of each temperature bin 
bysort avg_avg_temp_id: egen avg_bin=mean(avg_avg_temp)
table avg_bin, contents(sum trip_distance)

// Percentage change: ln form
gen ln_est_fe=ln(est_fe)
label variable ln_est_fe "ln(est_fe)"

xtreg ln_est_fe i(0/3 5/9)bn.avg_avg_temp_id, fe robust cluster(brand_group)
eststo m1
xtreg ln_est_fe i(0/3 5/9)bn.avg_avg_temp_id i.year_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo m2
xtreg ln_est_fe i(0/3 5/9)bn.avg_avg_temp_id i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo m3

estadd local Year_fixed_effects "NO":m1
estadd local Year_fixed_effects "YES":m2 m3
estadd local Month_fixed_effects "NO":m1 m2
estadd local Month_fixed_effects "YES":m3
estadd local Vehicle_fixed_effects "YES": m1 m2 m3

esttab m1 m2 m3 using ac_table1.csv, label b(3)  se(3) ar2 drop(*.month_mid *.year_mid) noconstant ///
		compress star(* 0.1 ** 0.05 *** 0.01) ///
		scalars(Year_fixed_effects Month_fixed_effects Vehicle_fixed_effects)  append   

		
/*==================================================================================*/
/*		          Robustness check 1: different temperature bins	        		*/
/*==================================================================================*/
// The second set of bins: 13 bins with 3˚C as an interval; the baseline group is 19-22˚C 
egen bin2 = cut(avg_avg_temp), at(-30,-2,1,4,7,10,13,16,19,22,25,28,31,40) icodes 

xtreg est_fe i(0/7 9/12)bn.bin2, fe robust cluster(brand_group)
eststo l1
xtreg est_fe i(0/7 9/12)bn.bin2 i.year_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo l2
xtreg est_fe i(0/7 9/12)bn.bin2 i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
eststo l3

estadd local Year_fixed_effects "NO":l1
estadd local Year_fixed_effects "YES":l2 l3
estadd local Month_fixed_effects "NO":l1 l2
estadd local Month_fixed_effects "YES":l3
estadd local Vehicle_fixed_effects "YES": l1 l2 l3

esttab l1 l2 l3 using ac_table1.csv, label b(3)  se(3) ar2 drop(*.month_mid *.year_mid) noconstant ///
		compress star(* 0.1 ** 0.05 *** 0.01) ///
		scalars(Year_fixed_effects Month_fixed_effects Vehicle_fixed_effects) ///
		append

// Get the weighting of each temperature bin 
bysort bin2: egen avg_bin2=mean(avg_avg_temp)
table avg_bin2, contents(sum trip_distance)


/*==================================================================================*/
/*		          Robustness check 2: cubic spline regression	         			*/
/*==================================================================================*/

// nknots=5
mkspline2  temp_spline1 = avg_avg_temp, cubic nknots(5)  displayknots
eststo m1: xtreg est_fe temp_spline1* i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
adjustrcspline, link(identity) ///
		yline(0) ciopts(color(eltblue*0.5)) lineopts(color(eltblue*1.5)) ///
		xtitle("Average temperature (˚C)") ytitle("Fuel consumption (L/100km)")  ///
		generate(y_spline1 lower upper) 

preserve
collapse (mean) avg_avg_temp (mean) y_spline1  (mean) lower (mean) upper, by(avg_avg_temp_id)   // mean predictions
egen y_min=min(y_spline1)
gen increase=y_spline1-y_min
gen lb=lower -y_min
gen ub=upper -y_min
export delimited using "spline1.csv", replace
restore

// nknots=7
mkspline2  temp_spline2 = avg_avg_temp, cubic nknots(7)  displayknots
eststo m2: xtreg est_fe temp_spline2* i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)
drop lower upper
adjustrcspline, link(identity) ///
		yline(0) ciopts(color(eltblue*0.5)) lineopts(color(eltblue*1.5)) ///
		xtitle("Average temperature (˚C)") ytitle("Fuel consumption (L/100km)")  ///
		generate(y_spline2 lower upper) 

preserve
collapse (mean) avg_avg_temp (mean) y_spline2  (mean) lower (mean) upper, by(avg_avg_temp_id)   // mean predictions
egen y_min=min(y_spline2)
gen increase=y_spline2-y_min
gen lb=lower -y_min
gen ub=upper -y_min
export delimited using "spline2.csv", replace
restore	
		
// Overall fuel consumption increase
* nknots=5
egen y_min=min(y_spline1)
gen fuel_increase=y_spline1-y_min
sort avg_avg_temp
bysort avg_avg_temp: gen unique=1 if _n==1  //  108,194	
preserve
collapse (mean) fuel_increase (sum) trip_distance, by(avg_avg_temp)
sort avg_avg_temp
egen total_distance=sum(trip_distance)
gen weighting=trip_distance/total_distance
egen results=sum(fuel_increase*weighting) if avg_avg_temp>=25
sum results, detail
restore

* nknots=7
drop y_min fuel_increase
egen y_min=min(y_spline2)
gen fuel_increase=y_spline2-y_min
sort avg_avg_temp
bysort avg_avg_temp: gen unique=1 if _n==1  //  108,194	
preserve
collapse (mean) fuel_increase (sum) trip_distance, by(avg_avg_temp)
sort avg_avg_temp
egen total_distance=sum(trip_distance)
gen weighting=trip_distance/total_distance
egen results=sum(fuel_increase*weighting) if avg_avg_temp>=25
sum results, detail
restore

/*==================================================================================*/
/*		                            Plot figures	                    			*/
/*==================================================================================*/
// "ac_figure1.csv" is a file compiling all previuos regression results
import delimited using ac_figure1.csv, clear

// 
graph twoway rarea upper lower avg_bin if form=="natural"&bin=="avg_avg_temp_id", color(eltblue*0.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="avg_avg_temp_id", color(eltblue*1.5) || ///
			 scatter mean_effect avg_bin if form=="natural"&bin=="avg_avg_temp_id", color(eltblue*1.5) msymbol(O) ///
			 xtitle("Average temperature (˚C)") ytitle("Percentage change of fuel consumption (%)") ///
			 legend(off) xlabel(-5 0 5 10 15 20 25 30 ) ylabel(0 "0.0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1.0 "1.0")
graph export bin1.png, replace

graph twoway rarea upper lower avg_bin if form=="ln"&bin=="avg_avg_temp_id", color(eltblue*0.5) || ///
			 line mean_effect avg_bin if form=="ln"&bin=="avg_avg_temp_id", color(eltblue*1.5) || ///
			 scatter mean_effect avg_bin if form=="ln"&bin=="avg_avg_temp_id", color(eltblue*1.5) msymbol(O) ///
			 xtitle("Average temperature (˚C)") ytitle("Percentage change of fuel consumption (%)") ///
			 legend(off) xlabel(-5 0 5 10 15 20 25 30 ) ylabel(0 "0" 0.02 "2" 0.04 "4" 0.06 "6" 0.08 "8" 0.10 "10")
graph export bin1_percentage.png, replace

graph twoway rarea upper lower avg_bin if form=="natural"&bin=="bin2", color(eltblue*0.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="bin2", color(eltblue*1.5)  || ///
			 scatter mean_effect avg_bin if form=="natural"&bin=="bin2", color(eltblue*1.5) msymbol(O) ///
			 xtitle("Average temperature (˚C)") ytitle("Fuel consumption increase (L/100km)") ///
			 legend(off) xlabel(-5 0 5 10 15 20 25 30 ) ylabel(0 "0.0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1.0 "1.0")
graph export bin2.png, replace

graph twoway rarea upper lower avg_bin if form=="natural"&bin=="avg_avg_temp_id", color(eltblue*0.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="avg_avg_temp_id", color(eltblue*1.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="bin2", color(dkorange) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="spline1", lpattern(dash) color(red)  || ///
			 line mean_effect avg_bin if form=="natural"&bin=="spline2", lpattern(dash) color(dkorange) ///
			 , legend(label(1 "95% CI") label(2 "bin1") label( 3 "bin2") label( 4 "cubic, 5 knots") label( 5 "cubic, 7 knots")) ///
			  xtitle("Average temperature (˚C)") ytitle("Fuel consumption increase (L/100km)")	 ///
			  xlabel(-5 0 5 10 15 20 25 30 ) ylabel(0 "0.0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1.0 "1.0")
graph export multiple_lines.png, replace

graph twoway rarea upper lower avg_bin if form=="natural"&bin=="avg_avg_temp_id"&avg_bin<=36&avg_bin>=15 , color(eltblue*0.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="avg_avg_temp_id"&avg_bin<=36&avg_bin>=15 , color(eltblue*1.5) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="bin2"&avg_bin<=36&avg_bin>=15 , color(dkorange) || ///
			 line mean_effect avg_bin if form=="natural"&bin=="spline2"&avg_bin<=36&avg_bin>=15 , lpattern(dash) color(red)  || ///
			 line mean_effect avg_bin if form=="natural"&bin=="spline3"&avg_bin<=36&avg_bin>=15 , lpattern(dash) color(dkorange) ///
			 , legend(label(1 "95% CI") label(2 "bin1") label( 3 "bin2") label( 4 "cubic, 5 knots") label( 5 "cubic, 7 knots")) ///
			  xtitle("Average temperature (˚C)") ytitle("Fuel consumption increase (L/100km)") ///
			  xscale(range(20 36))	 xlabel(15 20 25 30 35) ylabel(0 "0.0" 0.2 "0.2" 0.4 "0.4" 0.6 "0.6" 0.8 "0.8" 1.0 "1.0")
graph export multiple_lines_hot.png, replace


/*==================================================================================*/
/*		      Heterogenenous temperature responses by car models     			    */
/*==================================================================================*/
xtset uuid edate_hour

xtreg est_fe i.brand_group#i(1/4 6/10)bn.avg_avg_temp_id i.year_mid ib2.month_mid num_gap_day avg_humid avg_pcp, fe robust cluster(brand_group)

// Export results
matrix coeff = e(b)
matrix vc = e(V)
local df_r = e(df_r)
contract brand_group brand_displacement mean_est_fe

quiet gen coeff_brand = .
quiet gen pval_brand = .

local tmax = 9
local tmax_p = 10
forval t = 1/`tmax_p'{
	quiet gen coeff_brand_temp_`t' = .
	quiet gen pval_brand_temp_`t' = .
}

forval g = 1/`gmax' {
	forval t = 1/`tmax_p'{
		quiet replace coeff_brand_temp_`t' = coeff[1,(`g'-1)*`tmax'+`t'] if `t' < 5 & brand_group == `g'
		quiet replace pval_brand_temp_`t' = 2 * ttail(`df_r',abs(coeff[1,(`g'-1)*`tmax'+`t'] / sqrt(vc[(`g'-1)*`tmax'+`t',(`g'-1)*`tmax'+`t'])))  if `t' < 5 & brand_group == `g'
		quiet replace coeff_brand_temp_`t' = coeff[1,(`g'-1)*`tmax'+`t'-1] if `t' > 5 & brand_group == `g'
		quiet replace pval_brand_temp_`t' = 2 * ttail(`df_r',abs(coeff[1,(`g'-1)*`tmax'+`t'-1] / sqrt(vc[(`g'-1)*`tmax'+`t'-1,(`g'-1)*`tmax'+`t'-1])))  if `t' > 5 & brand_group == `g'
		quiet replace coeff_brand_temp_`t' = 0 if `t' == 5
		quiet replace pval_brand_temp_`t' = . if `t' == 5
	}
}

order brand_* coeff_* pval_*
outsheet brand_* mean_est_fe coeff_brand_temp* pval_brand_temp* using result_fe.csv, comma replace

