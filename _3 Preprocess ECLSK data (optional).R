# ----------------------------------------------------
# Methodology and Statistics in Psychology
# Master Thesis
# Maaike Jorink, s1907387
# ----------------------------------------------------

##### Load ECLSK data ######

library( data.table )
ICPSR <- fread( "28023-0001-Data.tsv" )

## Load theta errata data:
theta_err <- read.csv2( "eclsk_theta_errata.csv", header = T )

## Replace old scale values with revised theta values:
# Reading
old_names_read <- c( "C2R4RSCL", "C4R4RSCL", "C5R4RSCL", "C6R4RSCL", "C7R4RSCL" )
errata_names_read <- c( "C2R4RTHT_R", "C4R4RTHT_R", "C5R4RTHT_R", "C6R4RTHT_R", "C7R4RTHT_R" )
# Math
old_names_math <- c( "C2R4MTSC", "C4R4MTSC", "C5R4MTSC", "C6R4MTSC", "C7R4MTSC" )
errata_names_math <- c( "C2R4MTHT_R", "C4R4MTHT_R", "C5R4MTHT_R", "C6R4MTHT_R", "C7R4MTHT_R" )
# Science
old_names_scie <- c( "C5R2STHT", "C6R2STHT", "C7R2STHT" )
errata_names_scie <- c( "C5R2STHT_R", "C6R2STHT_R", "C7R2STHT_R" )

# Add to dataset
ICPSR[, old_names_read ] <- theta_err[, errata_names_read ]
ICPSR[, old_names_math ] <- theta_err[, errata_names_math ]
ICPSR[, old_names_scie ] <- theta_err[, errata_names_scie ]

## Make dataframe long format
library( tidyr )
data_wide <- ICPSR[, c( "CHILDID", "DOBMM", "DOBYY", 
                        "C2ASMTMM", "C4ASMTMM", "C5ASMTMM", "C6ASMTMM", "C7ASMTMM",
                        "GENDER", "RACE", "WKSESL", "C1GMOTOR", "C1FMOTOR", "T1INTERN", 
                        "T1EXTERN", "T1INTERP", "T1CONTRO", "P1FIRKDG" ) ]

## Change all missings to NA
data_wide[ data_wide == -1 ] <- NA
data_wide[ data_wide == -2 ] <- NA
data_wide[ data_wide == -7 ] <- NA
data_wide[ data_wide == -8 ] <- NA
data_wide[ data_wide == -9 ] <- NA

## Compute number of months passed the first measurement occasion
data_wide[, "AGEC2" ] <- ( data_wide$C2ASMTMM - data_wide$C2ASMTMM )
data_wide[, "AGEC4" ] <- ( data_wide$AGEC2 + 12 + ( data_wide$C4ASMTMM - data_wide$C2ASMTMM ) )
data_wide[, "AGEC5" ] <- ( data_wide$AGEC4 + 24 + ( data_wide$C5ASMTMM - data_wide$C4ASMTMM ) )
data_wide[, "AGEC6" ] <- ( data_wide$AGEC5 + 24 + ( data_wide$C6ASMTMM - data_wide$C5ASMTMM ) )
data_wide[, "AGEC7" ] <- ( data_wide$AGEC6 + 36 + ( data_wide$C7ASMTMM - data_wide$C6ASMTMM ) )
data_wide[, "AGEBASELINE" ] <- ( 1999 - data_wide$DOBYY ) * 12 + data_wide$C2ASMTMM - data_wide$DOBMM
# Seperate variables for science as science only has three measurement moments
data_wide[, "AGEC_5" ] <- ( data_wide$C5ASMTMM - data_wide$C5ASMTMM )
data_wide[, "AGEC_6" ] <- ( data_wide$AGEC_5 + 24 + ( data_wide$C6ASMTMM - data_wide$C5ASMTMM ) )
data_wide[, "AGEC_7" ] <- ( data_wide$AGEC_6 + 36 + ( data_wide$C7ASMTMM - data_wide$C6ASMTMM ) )

### Make seperate dataset per ability
## Reading 
data_wide_read <- data_wide
data_wide_read[, c( "C2R4RSCL", "C4R4RSCL", "C5R4RSCL", "C6R4RSCL", "C7R4RSCL" ) ] <-
  ICPSR[, c( "C2R4RSCL", "C4R4RSCL", "C5R4RSCL", "C6R4RSCL", "C7R4RSCL" ) ] # Add reading score variables
# Make data long form
data_long_read <- gather( data_wide_read, reading, score, C2R4RSCL:C7R4RSCL ) 
data_long_read2 <- gather( data_wide_read, asmtmm, months, AGEC2:AGEC7 ) 
orddata_read <- data_long_read[ order( data_long_read$CHILDID ), ]
orddata_read2 <- data_long_read2[ order( data_long_read2$CHILDID ), ]
orddata_read[, c( "asmtmm","months" ) ] <- orddata_read2[, c( "asmtmm", "months" ) ]

# datafile with only necessary variables.
newdata_read <- orddata_read[, c( "CHILDID", "GENDER", "RACE", "WKSESL", "C1GMOTOR", "C1FMOTOR",
                                  "T1INTERN", "T1EXTERN", "T1INTERP", "T1CONTRO", "P1FIRKDG", "AGEBASELINE", 
                                  "reading", "score", "asmtmm", "months" ) ]
# Use only cases that have a complete measurement moment
completedata_read <- newdata_read[ complete.cases( newdata_read ), ]

# Select IDs with at least 5 measurement occasions
select_ids_read <- names( table( completedata_read$CHILDID )[ table( completedata_read$CHILDID ) ==  5 ] )
readdata <- completedata_read[ completedata_read$CHILDID %in% select_ids_read, ]

# Make race an unordered factor
readdata$RACE <- factor( readdata$RACE, ordered = F )
# Transform months ^ (1/2)
readdata$months <- readdata$months ^ ( 1/2 )

## Math
data_wide_math <- data_wide
data_wide_math[, c( "C2R4MTSC", "C4R4MTSC", "C5R4MTSC", "C6R4MTSC", "C7R4MTSC" ) ] <-
  ICPSR[, c( "C2R4MTSC", "C4R4MTSC", "C5R4MTSC", "C6R4MTSC", "C7R4MTSC" ) ] # Add math score variables
# Make data long form
data_long_math <- gather( data_wide_math, math, score, C2R4MTSC:C7R4MTSC ) 
data_long_math2 <- gather( data_wide_math, asmtmm, months, AGEC2:AGEC7 ) 
orddata_math <- data_long_math[ order( data_long_math$CHILDID ), ]
orddata_math2 <- data_long_math2[ order( data_long_math2$CHILDID ), ]
orddata_math[, c( "asmtmm","months" ) ] <- orddata_math2[, c( "asmtmm", "months" ) ]

# datafile with only necessary variables.
newdata_math <- orddata_math[, c( "CHILDID", "GENDER", "RACE", "WKSESL", "C1GMOTOR", "C1FMOTOR",
                                  "T1INTERN", "T1EXTERN", "T1INTERP", "T1CONTRO", "P1FIRKDG", "AGEBASELINE", 
                                  "math", "score", "asmtmm", "months" ) ]
# Use only cases that have a complete measurement moment
completedata_math <- newdata_math[ complete.cases( newdata_math ), ]

# Select IDs with at least 5 measurement occasions
select_ids_math <- names( table( completedata_math$CHILDID )[ table( completedata_math$CHILDID ) ==  5 ] )
mathdata <- completedata_math[ completedata_math$CHILDID %in% select_ids_math, ]

# Make race an unordered factor
mathdata$RACE <- factor( mathdata$RACE, ordered = F )
# Transform months ^ (1/2)
mathdata$months <- mathdata$months ^ ( 1/2 )

## Science
data_wide_scie <- data_wide
data_wide_scie[, c( "C5R2STHT", "C6R2STHT", "C7R2STHT") ] <- 
  ICPSR[, c( "C5R2STHT", "C6R2STHT", "C7R2STHT") ] # Add science score variables
# Make data long form
data_long_scie <- gather( data_wide_scie, science, score, C5R2STHT:C7R2STHT ) 
data_long_scie2 <- gather( data_wide_scie, asmtmm, months, AGEC_5:AGEC_7 ) 
orddata_scie <- data_long_scie[ order( data_long_scie$CHILDID ), ]
orddata_scie2 <- data_long_scie2[ order( data_long_scie2$CHILDID ), ]
orddata_scie[, c( "asmtmm","months" ) ] <- orddata_scie2[, c( "asmtmm", "months" ) ]

# datafile with only necessary variables.
newdata_scie <- orddata_scie[, c( "CHILDID", "GENDER", "RACE", "WKSESL", "C1GMOTOR", "C1FMOTOR",
                                  "T1INTERN", "T1EXTERN", "T1INTERP", "T1CONTRO", "P1FIRKDG", "AGEBASELINE", 
                                  "science", "score", "asmtmm", "months" ) ]

# Use only cases that have a complete measurement moment
completedata_scie <- newdata_scie[ complete.cases( newdata_scie ), ]

# Select IDs with at least 3 measurement occasions
select_ids_scie <- names( table( completedata_scie$CHILDID )[ table( completedata_scie$CHILDID ) ==  3 ] )
sciedata <- completedata_scie[ completedata_scie$CHILDID %in% select_ids_scie, ]

# Make race an unordered factor
sciedata$RACE <- factor( sciedata$RACE, ordered = F )
# Transform months ^ (1/2)
sciedata$months <- sciedata$months ^ ( 2/3 )

## Save datasets
save(sciedata, file = "Science ability data.Rdata")
save(mathdata, file = "Math ability data.Rdata")
save(readdata, file = "Reading ability data.Rdata")