
********************************************************************************
*	Quantifying the contributions of active and past placental malaria infection 
*	to low birthweight in a high transmission area of Uganda
*	4_Participant Flowchart
********************************************************************************

clear all

* ------------------------------------
*	User inputs: 
* ------------------------------------

*** Set path directory below ***
capture cd "[insert path directory]"

* -------------------------------------
*	Number enrolled (n=2757)
* -------------------------------------
use "Databases/DPSP individual level analysis database_FINAL.dta", clear

keep id withdrawalreason 
merge 1:1 id using "Databases/DPSP delivery analysis database_FINAL.dta"

* -------------------------------------
*	Number delivered (n=2538)
* -------------------------------------
drop if _merge==1 // n=219

* -------------------------------------
*	Number of live births (n=2459)
* -------------------------------------
keep if infstatus==1 // n=79

* ------------------------------------------------------------
*	Number with placental histopathology results (n=2323)
* ------------------------------------------------------------
drop if rogerson==. // n=136

* ---------------------------------------------------
*	Number with birthweight data (n=2322)
* ---------------------------------------------------
drop if weight==. // n=1




