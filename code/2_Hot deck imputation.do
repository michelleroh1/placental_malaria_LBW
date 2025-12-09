
********************************************************************************
*	Quantifying the contributions of active and past placental malaria infection 
*	to low birthweight in a high transmission area of Uganda
*	2_Hot deck imputation 
********************************************************************************

clear all 

* ------------------------------------
*	User inputs: 
* ------------------------------------

*** Set path directory below ***
capture cd "[insert path directory]"

*** Install hotdeckvar package (if not already done so) ***
ssc install hotdeckvar

* ------------------------------------
*	Upload dataset for imputation
* ------------------------------------
use "Databases/Created Databases/DPSP Study - Delivery and Histopath Results.dta", clear

keep id activeHP propfibrincat anyplcmal ///
		   preterm SGA TLBW LBWdich ///
		   graviddich gravidcat ///
		   treatmentarm age_enroll age2 wealthcat bmi_cat bmi_enroll ///
		   GAenroll gaenroll2 educ_cat infstatus weight plcmal_cat male BSenroll SDdosecat
		   
* ---------------------------------------
*	Identify # of missing in each table 
* ---------------------------------------
		   
// Restrict data to analytic sample 
keep if infstatus == 1
keep if !missing(plcmal_cat)
keep if weight!=.
		   
// Variables included in adjusted models
local vars activeHP propfibrincat anyplcmal ///
		   preterm SGA LBWdich graviddich ///
		   treatmentarm age_enroll wealthcat bmi_cat bmi_enroll ///
		   GAenroll educ_cat male
		   
// Calculate % (n/N) missing for each variable 	
foreach var of local vars {
    display "----------------------------------------"
    quietly count if missing(`var')
    local miss = r(N)
    local total = _N
    local percent = 100 * `miss' / `total'
    display as text "`var': " ///
	 as result %4.1f `percent' ///
     as text "% missing (" `miss' "/" `total' ")"
    }

* ----------------------------------------------------
*	Use hot deck imputation given the small % of missingness (<0.4%) of variables 
* ----------------------------------------------------

set seed 123
hotdeckvar wealthcat bmi_enroll bmi_cat, suffix("_hd")

// imputed wealthcat variable
tab wealthcat wealthcat_hd, missing
tab bmi_cat bmi_cat_hd, missing

// Using imputed BMI variable, create BMI^2 variable: bmi_enroll_hd2
capture drop bmi_enroll_hd2
gen bmi_enroll_hd2 = bmi_enroll_hd^2

save "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed.dta", replace
		   