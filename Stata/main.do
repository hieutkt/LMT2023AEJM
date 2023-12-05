cd "C:\Users\Cristian\Documents\Replication Capital Allocation\PWT"

// Uploding Penn World Table version 10.01
use "pwt1001.dta", clear  

// Keep only peruvian data
keep if country == "Peru"

drop countrycode
save "temp.dta",replace 



// this data set contains some information about trade from the IMF 
use "peru_data.dta", clear
destring year, replace 
cap gen country = "Peru"


// Merge both
merge 1:1 country year using blah.dta

drop _merge
save "temp.dta",replace 




// Open the data set that contains information about imports
import excel "Peru's Imports from China.xlsx", sheet("Sheet1") firstrow clear
destring year, replace

//merge
merge 1:1 year using blah2.dta



// Mantain observation only between 2000 and 2023
drop if year<=2000 | year >=2023
drop country* _merge


// re-scale imports
replace imports = imports / 1000

// graph for imports
twoway (connected imports year if year >= 2000) , ///
title("Peru's imports from China'") xline(2020, lstyle(foreground)) /// 
scheme(s2mono) ytitle("Thousands of USD")
graph export imports.pdf, replace  
graph export imports.png, width(600) height(450) replace


// graph for annual average number of hours of people engaged in production
twoway (connected avh year if year >= 2000) ///
(lfit avh  year if year >= 2000, leg(off)),  ///
title("Average annual hours worked by persons engaged") /// 
scheme(s2mono) ytitle(" ")



// graph: Commodity terms of trade
twoway (connected CommodityTermofTrade year if year >= 2000) ///
,  xline(2020, lstyle(foreground)) ///
title("ToT over time") scheme(s2mono) ytitle("Terms of trade")

graph export ToT.pdf, replace 
graph export ToT.png, width(600) height(450) replace

// graph: Import price

twoway (connected pl_m  year if year >= 2000,  leg(off)) ///
(lfit pl_m  year if year >= 2000, leg(off)) , ///
title("Import price over time") scheme(s2mono) ytitle("Price level of imports") 

graph export import_price.pdf, replace
graph export import_price.png, width(600) height(450) replace



/*-----------------------------------------------------------
Open dataset that contains information about entry and 
exit rates on the whole economy and manufacture sector. Source: INEI, demografía empresarial
----------------------------------------------------------*/

// Whole economy 
import excel "Entry exit rates in Peru.xlsx", sheet("Whole Economy") firstrow clear
drop J K L

rename exitrate exitrate_w
rename entryrate entryrate_w
rename Stockfirms  stock_w
rename entry entry_w
rename exit exit_w


save "temp.dta", replace



// Open sheet that contains information of the manufacture sector
import excel "Entry exit rates in Peru.xlsx", sheet("Manufacture Sector") firstrow clear

drop K L

rename exitrate exitrate_m

rename entryrate entryrate_m
rename stockfirms stock_m
rename newfirms entry_m
rename exitfirms  exit_m

merge 1:1 date using potato.dta


// time variable
gen yqdate = yq(year, quarter)
format yqdate %tq
tsset yqdate


// creating yoy and period by period growth rate for the stock of 
// new firms, firms that exit, and the total stock
local varGrowth "entry exit stock"
foreach jj of local varGrowth{
    cap drop growth_rate_m4_`jj'
	cap drop growth_rate_w4_`jj'
    gen growth_rate_m4_`jj' = (`jj'_m/ L4.`jj'_m - 1)*100
	gen growth_rate_w4_`jj' = (`jj'_w/ L4.`jj'_w - 1)*100
	
	
	cap drop growth_rate_m1_`jj'
	cap drop growth_rate_w1_`jj'
    gen growth_rate_m1_`jj' = (`jj'_m/ L1.`jj'_m - 1)*100
	gen growth_rate_w1_`jj' = (`jj'_w/ L1.`jj'_w - 1)*100

}
 



// graph of the growth rate of the stock, year-to-year
twoway (connected growth_rate_w4_stock date if yqdate>=tq(2013q3)) ///
(connected growth_rate_m4_stock date if yqdate>=tq(2013q3)), /// 
xline(2020, lstyle(foreground)) ///
title("Growth rate of the stock of firms, year to year") scheme(s2mono) ytitle("Growth rate") ///
legend(order (1 "Whole Economy" 2 "Manufacture"))
graph export gw_yoy.pdf, replace
graph export gw_yoy.png, width(600) height(450) replace

// graph of the growth rate of the stock, period-to-period
twoway (connected growth_rate_w1_stock date if yqdate>=tq(2013q3)) ///
(connected growth_rate_m1_stock date if yqdate>=tq(2013q3)), /// 
xline(2020, lstyle(foreground)) ///
title("Growth rate of the stock of firms, next period") scheme(s2mono) ytitle("Growth rate") ///
legend(order (1 "Whole Economy" 2 "Manufacture"))
graph export gw_period.pdf, replace
graph export gw_period.png, width(600) height(450) replace

// graph: entry rate
twoway (connected entryrate_w date if yqdate>=tq(2013q3)) ///
(connected entryrate_m date if yqdate>=tq(2013q3)), /// 
xline(2020, lstyle(foreground)) ///
title("Whole economy vs Manufacture sector") scheme(s2mono) ytitle("Entry rate") ///
legend(order (1 "Whole Economy" 2 "Manufacture"))
graph export entry_rate.pdf, replace
graph export entry_rate.png, width(600) height(450) replace

// graph: exit rate
twoway (connected exitrate_w date if yqdate>=tq(2013q3)) ///
(connected exitrate_m date if yqdate>=tq(2013q3)), /// 
xline(2020, lstyle(foreground)) ///
title("Whole economy vs Manufacture sector") scheme(s2mono) ytitle("Exit rate") ///
legend(order (1 "Whole Economy" 2 "Manufacture"))
graph export exit_rate.pdf, replace
graph export exit_rate.png, width(600) height(450) replace
 







