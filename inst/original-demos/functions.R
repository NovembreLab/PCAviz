# Return PC1 after a rotation and scaling, followed by a translation.
# Inputs sx and sy are scaling factors; inputs x and y translate the
# co-ordinates; and input argument p specifies the rotation angle in
# radians.
rotPC1 <- function (p,sx = 1,sy = 1,x = 0,y = 0)
  cos(p)*(sy*PCA[,"PC1"] - y) - sin(p)*(sx*PCA[,"PC2"] - x)      

# Return PC2 after a rotation and scaling, followed by a translation.
# Inputs sx and sy are scaling factors; inputs x and y translate the
# co-ordinates; and input argument p specifies the rotation angle in
# radians.
rotPC2 <- function (p,sx = 1,sy = 1,x = 0,y = 0)
  cos(p)*(sx*PCA[,"PC2"] - x) + sin(p)*(sy*PCA[,"PC1"] - y)

# Fix up some of the country-of-origin annotations for the plots.
fix.country.names <- function (x) {
  x[grep("Scotland",x)]="United Kingdom"
  x[grep("Swiss-German",x)]="Switzerland"
  x[grep("Swiss-Italian",x)]="Switzerland"
  x[grep("Swiss-French",x)]="Switzerland"
  x[grep("Kosovo",x)]="Serbia and Montenegro"
  x[grep("Yugoslavia",x)]="Serbia and Montenegro"
  x[grep("Serbia",x)]="Serbia and Montenegro"
  x[grep("Bosnia",x)]="Bosnia and Herzegovina"
  x[grep("Bosnia-Herzegovina",x)]="Bosnia and Herzegovina"
  x[grep("Russia",x)]="Russian Federation"
  x[grep("Macedonia",x)]="Macedonia, The Former Yugoslav Republic Of"
  return(x)
}

# An alternative correction of the country-of-origin annotations that
# is better suited for drawing maps.
fix.country.names.maps <- function (x) {
  x[grep("Wales",x)]="United Kingdom"
  x[grep("Scotland",x)]="United Kingdom"
  x[grep("Yugoslavia",x)]="Serbia and Montenegro"
  x[grep("Serbia",x)]="Serbia and Montenegro"
  x[grep("Bosnia",x)]="Bosnia and Herzegovina"
  x[grep("Bosnia-Herzegovina",x)]="Bosnia and Herzegovina"
  x[grep("Kosovo",x)]="Serbia and Montenegro"
  x[grep("Russia",x)]="Russian Federation"
  x[grep("Macedonia",x)]="Macedonia, The Former Yugoslav Republic Of"
  x[grep("Swiss-German",x)]="Switzerland"
  x[grep("Swiss-Italian",x)]="Switzerland"
  x[grep("Swiss-French",x)]="Switzerland"
  return(x)
}

# An updated function for fixing up some of the country-of-origin
# annotations for the plots.
fix.country.names2 <- function (x) {
  x[grep("Scotland",x)]="United Kingdom"
  x[grep("Swiss-German",x)]="Switzerland"
  x[grep("Swiss-Italian",x)]="Switzerland"
  x[grep("Swiss-French",x)]="Switzerland"
  x[grep("Yugoslavia",x)]="Serbia"
  x[grep("Serbia and Montenegro",x)]="Serbia"
  x[grep("Kosovo",x)]="Serbia"
  x[grep("Bosnia",x)]="Bosnia and Herz."
  x[grep("Bosnia-Herzegovina",x)]="Bosnia and Herz."
  x[grep("Russian Federation",x)]="Russia"
  x[grep("Czech Republic",x)]="Czech Rep."
  x[grep("Macedonia, The Former Yugoslav Republic Of",x)]="Macedonia"
  return(x)
}

fix.country.abbr <- function (x) {
  x[grep("UK",x)]="GB"
  x[grep("YG",x)]="RS"
  x[grep("Sct",x)]="GB"
  x[grep("KS",x)]="RS"
  return(x)
}

abbrvlabels <- function (x) {
  x[grep("United Kingdom",x)]="GB"
  x[grep("Wales",x)]="Wales"
  x[grep("Scotland",x)]="Sct"
  x[grep("Czech Republic",x)]="CZ"
  x[grep("Yugoslavia",x)]="YG"
  x[grep("Bosnia and Herzegovina",x)]="BA"
  x[grep("Kosovo",x)]="KS"
  x[grep("Russia",x)]="RU"
  x[grep("Iran",x)]="IR"
  x[grep("Macedonia",x)]="MK"
  x[grep("Swiss-German",x)]="CH"
  x[grep("Swiss-Italian",x)]="CH"
  x[grep("Swiss-French",x)]="CH"
  x[grep("Switzerland",x)]="CH"
  x[grep("Spain",x)]="ES"
  x[grep("Italy",x)]="IT"
  x[grep("Portugal",x)]="PT"
  x[grep("France",x)]="FR"
  x[grep("Germany",x)]="DE"
  x[grep("Greece",x)]="GR"
  x[grep("Ireland",x)]="IE"
  x[grep("Romania",x)]="RO"
  x[grep("Austria",x)]="AT"
  x[grep("Hungary",x)]="HU"
  x[grep("Belgium",x)]="BE"
  x[grep("Poland",x)]="PL"
  x[grep("Netherlands",x)]="NL"
  x[grep("Norway",x)]="NO"
  x[grep("Sweden",x)]="SE"
  x[grep("Finland",x)]="FI"
  x[grep("Latvia",x)]="LV"
  x[grep("Turkey",x)]="TR"
  x[grep("Croatia",x)]="HR"
  x[grep("Albania",x)]="AL"
  x[grep("Ukraine",x)]="UA"
  x[grep("Bulgaria",x)]="BG"
  x[grep("Slovenia",x)]="SI"
  x[grep("Slovakia",x)]="SK"
  x[grep("Cyprus",x)]="CY"
  x[grep("Denmark",x)]="DK"
  x[grep("Serbia",x)]="RS"
  return(x)
}

abbrvlabels2 <- function (x) {
  x[grep("United Kingdom",x)]="GB"
  x[grep("Wales",x)]="GB"
  x[grep("Scotland",x)]="GB"
  x[grep("Czech Republic",x)]="CZ"
  x[grep("Yugoslavia",x)]="RS"
  x[grep("Bosnia and Herzegovina",x)]="BA"
  x[grep("Kosovo",x)]="RS"
  x[grep("Russia",x)]="RU"
  x[grep("Iran",x)]="IR"
  x[grep("Macedonia",x)]="MK"
  x[grep("Swiss-German",x)]="CH"
  x[grep("Swiss-Italian",x)]="CH"
  x[grep("Swiss-French",x)]="CH"
  x[grep("Switzerland",x)]="CH"
  x[grep("Spain",x)]="ES"
  x[grep("Italy",x)]="IT"
  x[grep("Portugal",x)]="PT"
  x[grep("France",x)]="FR"
  x[grep("Germany",x)]="DE"
  x[grep("Greece",x)]="GR"
  x[grep("Ireland",x)]="IE"
  x[grep("Romania",x)]="RO"
  x[grep("Austria",x)]="AT"
  x[grep("Hungary",x)]="HU"
  x[grep("Belgium",x)]="BE"
  x[grep("Poland",x)]="PL"
  x[grep("Netherlands",x)]="NL"
  x[grep("Norway",x)]="NO"
  x[grep("Sweden",x)]="SE"
  x[grep("Finland",x)]="FI"
  x[grep("Latvia",x)]="LV"
  x[grep("Turkey",x)]="TR"
  x[grep("Croatia",x)]="HR"
  x[grep("Albania",x)]="AL"
  x[grep("Ukraine",x)]="UA"
  x[grep("Bulgaria",x)]="BG"
  x[grep("Slovenia",x)]="SI"
  x[grep("Slovakia",x)]="SK"
  x[grep("Cyprus",x)]="CY"
  x[grep("Denmark",x)]="DK"
  x[grep("Serbia",x)]="RS"
  return(x)
}

# Convert full population labels to abbreviated labels for plotting
# the Swiss samples.
abbrvlabels.swiss = function (x) {
  x[grep("United Kingdom",x)]="GB"
  x[grep("Wales",x)]="Wales"
  x[grep("Scotland",x)]="Sct"
  x[grep("Czech Republic",x)]="CZ"
  x[grep("Yugoslavia",x)]="YG"
  x[grep("Bosnia and Herzegovina",x)]="BA"
  x[grep("Kosovo",x)]="KS"
  x[grep("Russia",x)]="RU"
  x[grep("Iran",x)]="IR"
  x[grep("Macedonia",x)]="MK"
  x[grep("Swiss-German",x)]="CH-G"
  x[grep("Swiss-Italian",x)]="CH-I"
  x[grep("Swiss-French",x)]="CH-F"
  x[grep("Switzerland",x)]="CH"
  x[grep("Spain",x)]="ES"
  x[grep("Italy",x)]="IT"
  x[grep("Portugal",x)]="PT"
  x[grep("France",x)]="FR"
  x[grep("Germany",x)]="DE"
  x[grep("Greece",x)]="GR"
  x[grep("Ireland",x)]="IE"
  x[grep("Romania",x)]="RO"
  x[grep("Austria",x)]="AT"
  x[grep("Hungary",x)]="HU"
  x[grep("Belgium",x)]="BE"
  x[grep("Poland",x)]="PL"
  x[grep("Netherlands",x)]="NL"
  x[grep("Norway",x)]="NO"
  x[grep("Sweden",x)]="SE"
  x[grep("Finland",x)]="FI"
  x[grep("Latvia",x)]="LV"
  x[grep("Turkey",x)]="TR"
  x[grep("Croatia",x)]="HR"
  x[grep("Albania",x)]="AL"
  x[grep("Ukraine",x)]="UA"
  x[grep("Bulgaria",x)]="BG"
  x[grep("Slovenia",x)]="SI"
  x[grep("Slovakia",x)]="SK"
  x[grep("Cyprus",x)]="CY"
  x[grep("Denmark",x)]="DK"
  x[grep("Serbia",x)]="RS"
  return(x)
}
