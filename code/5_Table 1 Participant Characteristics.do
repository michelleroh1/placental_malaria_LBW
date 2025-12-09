
*********************************************************************************
*	Quantifying the contributions of active and past placental malaria infection 
* to low birthweight in a high transmission area of Uganda
*	5_Table 1 Participant Characteristics 
*********************************************************************************

clear all

* ------------------------------------
*	User inputs: 
* ------------------------------------

*** Set path directory below ***
capture cd "[insert path directory]"

* ------------------------------------
*	Upload dataset
* ------------------------------------
use "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta", clear

* ----------------------------------------------------
*	Ensure data is restricted to analytic study sample
* ----------------------------------------------------
keep if infstatus == 1
keep if !missing(plcmal_cat)
keep if weight!=.


* -----------------------------------------------------------------------------
*	Generate Table 1 for Active Infection and Past Infection Separately 
* -----------------------------------------------------------------------------

/*	Active Infection --------- */

table1_mc, ///
	by(activeHP) ///
	vars( ///
		 treatmentarm cat %4.1f \ ///
		 graviddich cat %4.1f \ ///
		 GAenroll contn %4.1f \ ///
		 SDdosecat cat %4.1f \ ///
		 age_enroll contn %4.1f \ ///
		 bmi_cat cat %4.1f \ ///
		 wealthcat cat %4.1f \ ///
		 educ_cat cat %4.1f ) ///
	nospace onecol total(after) ///
	saving("Results/Table 1 - Active Infection Exposure Levels.xlsx", replace)
	
/*	Past Infection --------- */

table1_mc, ///
	by(propfibrincat) ///
	vars( ///
		 treatmentarm cat %4.1f \ ///
		 graviddich cat %4.1f \ ///
		 GAenroll contn %4.1f \ ///
		 SDdosecat cat %4.1f \ ///
		 age_enroll contn %4.1f \ ///
		 bmi_cat cat %4.1f \ ///
		 wealthcat cat %4.1f \ ///
		 educ_cat cat %4.1f ) ///
	nospace onecol total(after) ///
	saving("Results/Table 1 - Pigment Exposure Levels.xlsx", replace)


	
	
