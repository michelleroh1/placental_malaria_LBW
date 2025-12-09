
**********************************************************************
*	Quantifying the contributions of active and past placental malaria infection 
*	to low birthweight in a high transmission area of Uganda
*	1_Clean dataset 
**********************************************************************

clear all 

* ------------------------------------
*	User inputs: 
* ------------------------------------

*** Set path directory below ***
* Set the path to the main project directory
local projectdir "[insert path to project directory]"
capture cd "`projectdir'"

*** Create required subfolders if they do not already exist ***
capture mkdir "`projectdir'/Databases" // insert raw data files here
capture mkdir "`projectdir'/Databases/Created Databases" 

* --------------------------------------------------------
*	Merge delivery + enrollment database		
* --------------------------------------------------------

use "Databases/DPSP delivery analysis database_FINAL.dta", clear // Note files are located in a "Databases" folder within the path directory 

tempfile delivery
save "`delivery'"

use "Databases/DPSP enrollment analysis database_FINAL.dta", clear 

keep id GAenroll netlastnight weight height educ wealthcat AGE GAenroll BSdich
rename (weight BSdich) (enroll_wt BSenroll)

merge 1:1 id using "`delivery'"
keep if _merge==3
drop _merge

* --------------------------------------------------------
*	Create variables 		
* --------------------------------------------------------

// 3-level education variable
capture drop educ_cat
gen educ_cat = educ
replace educ_cat = 2 if educ>=3
replace educ_cat = . if educ==.
label define educ_cat_label 0 "None" 1 "Primary" 2 "O level and above"
label values educ_cat educ_cat_label

// Binary variable for primigravidae 
gen primigravid = (gravidcat ==1)
replace primigravid = . if gravidcat==.

// Infant sex male (Y/N)? 
capture drop male 
gen male = (gender==1)
replace male = . if gender == . 

rename (AGE enroll_wt height) (age_enroll weight_enroll height_enroll)

// Age 
capture drop age2
gen age2 = age_enroll^2

// Gestational age at enrollment 
capture drop gaenroll2
gen gaenroll2 = GAenroll^2

// Generate BMI from weight and height at enrollment 
capture drop bmi_enroll 
gen bmi_enroll = weight_enroll/((height_enroll/100)^2)

// 4-level BMI at enrollment variable 
capture drop bmi_cat
gen bmi_cat = 0 if bmi_enroll < 18.5
replace bmi_cat = 1 if bmi_enroll >=18.5 & bmi_enroll <25
replace bmi_cat = 2 if bmi_enroll >=25 & bmi_enroll <30
replace bmi_cat = 3 if bmi_enroll >=30
replace bmi_cat = . if bmi_enroll == . 

label define cat_bmi 0 "<18.5" 1 ">=18.5% to <25" 2 ">=25 to <30" 3 ">=30"
label values bmi_cat cat_bmi

// Gestational week when study drugs were initiated 
capture drop SDdosecat 
gen SDdosecat = (GA1stdoseSD >=19)
replace SDdosecat = 2 if GA1stdoseSD == 24
label variable SDdosecat "Gestational weeks when study drugs started"
label define label_SDdosecat 0 "~16 GW" 1 "~20 GW" 2 "24 GW" 3 "No IPTp given"
label value SDdosecat label_SDdosecat
tab GA1stdoseSD SDdosecat, missing

// Categorize past infection into 3-level variable: none, mild (<10%), moderate (10 to <30%), severe (>=30%)
capture drop propfibrincat
gen propfibrincat = 0 if rogerson==0
replace propfibrincat = 1 if propfibrin >0 & propfibrin <.10
replace propfibrincat = 2 if propfibrin >=.10 & propfibrin <.30 
replace propfibrincat = 3 if propfibrin >=.30 & propfibrin!=.
replace propfibrincat = . if propfibrin ==. & rogerson!=0
replace propfibrincat = 0 if rogerson == 1 

label define fibrincat 0 "None" 1 "<10%" 2 "10 to <30%" 3 ">=30%"
label values propfibrincat fibrincat 

// 4-level categorical variable for placental malaria: 0 none, 1 past inf with mild pigment, 2 past inf with moderate pigment, 3 past inf with severe pigment, 4 any evidence of parasites (w/ or w/o pigment); values 1-3 do not include women with active infeciton.
capture drop plcmal_cat
gen plcmal_cat = 0 if rogerson == 0 
replace plcmal_cat = 1 if rogerson == 4 & propfibrincat==1
replace plcmal_cat = 2 if rogerson == 4 & propfibrincat==2
replace plcmal_cat = 3 if rogerson == 4 & propfibrincat==3
replace plcmal_cat = 4 if rogerson == 1 | rogerson == 2 | rogerson == 3
replace plcmal_cat = . if rogerson == . 

label define plcmal_label 0 "No Pigment" 1 "Mild Pigment" 2 "Moderate Pigment" 3 "Severe Pigment" 4 "Any Parasites (w or w/o Pigment)"
label value plcmal_cat plcmal_label

// Binary variable describing any parasites or pigment
capture drop anyplcmal
gen anyplcmal = (propfibrincat!=0 | activeHP==1)
replace anyplcmal = . if propfibrincat==. & activeHP==. 

// Term LBW (LBW among term infants; set preterm as missing)
capture drop TLBW
gen TLBW = (LBWdich == 1 & preterm !=1)
replace TLBW = . if infstatus!=1
replace TLBW = . if preterm==1

save "Databases/Created Databases/DPSP Study - Delivery and Histopath Results.dta", replace










