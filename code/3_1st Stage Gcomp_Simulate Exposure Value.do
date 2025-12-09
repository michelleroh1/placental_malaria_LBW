
********************************************************************************
*	Quantifying the contributions of active and past placental malaria infection 
*	to low birthweight in a high transmission area of Uganda
*	3_1st Stage Gcomp_Simulate Exposure Value
********************************************************************************

clear all 

* --------------------------------------
*	User inputs 
* --------------------------------------

*** Set path directory below ***
capture cd "[insert path directory]"


* --------------------------------------------------
*	Generate simulated exposure value of active infection 
* --------------------------------------------------

/*	1. Load dataset ------------------ */
	use "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed.dta", clear
	set seed 123

/*	2. Generate exposure distribution under counterfactual IPTp-SP scenario ------ */

	// These are conditional probabilities based on covariate distribtion, except
	// we fix IPTp to SP for all participants 

	// First create an identical dataset: 
	//		1 = original (the dataset used to build model)
	//		2 = exposed at the predicted level under IPTp-SP */ 
	expand 2, gen(cf_num)

	label define label_cf_num 0 "original dataset" 1 "all SP, sim activeHP"
	label values cf_num label_cf_num

	// In cf_num 2  databases, create the hypothetical scenario of IPTp-SP
		// (1) Drop observed activeHP 
		replace activeHP = . if cf_num == 1 
		// (2) Assume everyone would have received IPTp-SP
		replace treatmentarm = 1 if cf_num == 1
	
/*	3. Generate 'exposure' model to estimate counterfactual activeHP values ---------------------- */

	logistic activeHP i.treatmentarm##i.graviddich c.age_enroll age2 i.wealthcat_hd ///
					  c.bmi_enroll_hd c.bmi_enroll_hd2 GAenroll gaenroll2 i.educ_cat i.male if cf_num==0
	predict activeHP_cf if cf_num==1, pr			

/*	4. Simulate each woman's activeHP value based on predicted probability  --------- */
	gen activeHP_sim = rbinomial(1, activeHP_cf) if cf_num == 1

/*	5. Compare activeHP prevalence in simulated dataset to original dataset (but only in the IPTp-SP arm) --------- */
	
	// *** OVERALL *** // 
	tabstat activeHP activeHP_sim if treatmentarm==1, by(cf_num) statistic(mean)

	// *** PRIMIGRAVIDAE *** //
	tabstat activeHP activeHP_sim if graviddich==1 & treatmentarm==1, by(cf_num) statistic(mean)

	// *** MULTIGRAVIDAE *** // 
	tabstat activeHP activeHP_sim if graviddich==2 & treatmentarm==1, by(cf_num) statistic(mean)
	
/*	6. Keep counterfactual dataset only ------------ */
	keep if cf_num==1
	keep id activeHP_sim 

	merge 1:1 id using "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed.dta"
	drop _merge

	save "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp.dta", replace	
	
	
	
* --------------------------------------------------
*	Generate simulated exposure value for pigmentation 
* --------------------------------------------------
	
/*	1. Load dataset ------------------ */
	use "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed.dta", clear
	set seed 123

/*	2. Duplicate dataset: original (cf_num=0) and counterfactual (cf_num=1) ---------- */	
	expand 2, gen(cf_num)

	label define label_cf_num 0 "original dataset" 1 "all SP, sim propfibrincat"
	label values cf_num label_cf_num

	// In cf_num==1, drop observed propfibrincat and set treatment to IPTp-SP //
	replace propfibrincat = . if cf_num == 1
	replace treatmentarm = 1 if cf_num == 1

/*	3. Generate 'exposure model' (multinomial logistic regression) to estimate propfibrincat ---------- */
	mlogit propfibrincat i.treatmentarm##i.graviddich c.age_enroll age2 ///
		   i.wealthcat_hd c.bmi_enroll_hd c.bmi_enroll_hd2 GAenroll gaenroll2 i.educ_cat i.male if cf_num==0

/*	4. Get predicted probabilities for each category in the counterfactual dataset */
	predict p1 p2 p3 p4 if cf_num==1, pr
		   
/*	5. Simulate each woman's category based on probabilities ----------- */
	gen u = runiform() if cf_num==1

	gen propfibrincat_sim = .
	replace propfibrincat_sim = 0 if cf_num==1 & u <= p1
	replace propfibrincat_sim = 1 if cf_num==1 & u > p1 & u <= p1+p2
	replace propfibrincat_sim = 2 if cf_num==1 & u > p1+p2 & u <= p1+p2+p3
	replace propfibrincat_sim = 3 if cf_num==1 & u > p1+p2+p3

/*	6. Keep counterfactual dataset only ----------- */
	keep if cf_num==1
	keep id propfibrincat_sim 

	merge 1:1 id using "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp.dta"
	drop _merge

	save "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp.dta", replace

	

* -----------------------------------------------------------------------
*	Generate simulated exposure value for any parasites or pigment
* -----------------------------------------------------------------------

/*	1. Load dataset ------------------ */
	use "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp.dta", clear
	set seed 123

/*	2. Duplicate dataset: original (cf_num=0) and counterfactual (cf_num=1) ------ */
	expand 2, gen(cf_num)

	label define label_cf_num 0 "original dataset" 1 "all SP, sim anyplcmal"
	label values cf_num label_cf_num

	// In cf_num 2  databases, create the hypothetical scenario of IPTp-SP
		// (1) Drop observed activeHP 
		replace anyplcmal = . if cf_num == 1 
		// (2) Assume everyone would have received IPTp-SP
		replace treatmentarm = 1 if cf_num == 1

/*	3. Generate 'exposure' model to estimate counterfactual anyplcmal values ------------ */
	logistic anyplcmal i.treatmentarm##i.graviddich c.age_enroll age2 i.wealthcat_hd ///
					  c.bmi_enroll_hd c.bmi_enroll_hd2 GAenroll gaenroll2 i.educ_cat i.male if cf_num==0
	predict anyplcmal_cf if cf_num==1, pr			

/*	4. Simulate each woman's anyplcmal value based on predicted probability  --------- */
	gen anyplcmal_sim = rbinomial(1, anyplcmal_cf) if cf_num == 1

/*	5. Compare activeHP prevalence in simulated dataset to original dataset (but only in the IPTp-SP arm) --------- */
	
	// *** OVERALL *** // 
	tabstat anyplcmal anyplcmal_sim if treatmentarm==1, by(cf_num) statistic(mean)

	// *** PRIMIGRAVIDAE *** //
	tabstat anyplcmal anyplcmal_sim if graviddich==1 & treatmentarm==1, by(cf_num) statistic(mean)

	// *** MULTIGRAVIDAE *** // 
	tabstat anyplcmal anyplcmal_sim if graviddich==2 & treatmentarm==1, by(cf_num) statistic(mean)

/*	6. Keep counterfactual dataset only ------------ */
	keep if cf_num==1
	keep id anyplcmal_sim 

	merge 1:1 id using "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp.dta"
	drop _merge

	save "Databases/Created Databases/DPSP Study - Delivery and Histopath Results_HDimputed_Sim_Gcomp_FINAL.dta", replace
