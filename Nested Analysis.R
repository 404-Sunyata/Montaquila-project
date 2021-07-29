###Montaquila et al. 2021 Community Assembly of Earthworms###


### Code for Anlaysing Nestedness Ordered by Decreasing Fill, 
#Spearman Correlation with variables, and Fox's rule


# Adding vegan to library with nestednodf function
library(vegan)
library(Rfast)

`nestednodf` <- 
  function(comm, order = TRUE, weighted = FALSE, wbinary = FALSE) 
  {
    bin.comm <- ifelse(comm > 0, 1, 0)
    rfill <- rowSums(bin.comm)
    cfill <- colSums(bin.comm)
    if (!weighted)
      comm <- bin.comm
    if (order) {
      if (weighted) {
        rgrad <- rowSums(comm)
        cgrad <- colSums(comm)
        rorder <- order(rfill, rgrad, decreasing = TRUE)
        corder <- order(cfill, cgrad, decreasing = TRUE)
      } else {
        rorder <- order(rfill, decreasing = TRUE)
        corder <- order(cfill, decreasing = TRUE)
      }
      comm <- comm[rorder, corder]
      rfill <- rfill[rorder]
      cfill <- cfill[corder]
    }
    nr <- NROW(comm)
    nc <- NCOL(comm)
    fill <- sum(rfill)/prod(dim(comm))
    N.paired.rows <- numeric(nr * (nr - 1)/2)
    N.paired.cols <- numeric(nc * (nc - 1)/2)
    counter <- 0
    for (i in 1:(nr - 1)) {
      first <- comm[i, ]
      for (j in (i + 1):nr) {
        counter <- counter + 1
        if (rfill[i] <= rfill[j] || any(rfill[c(i, j)] == 0)) 
          next
        if (weighted) {
          second <- comm[j, ]
          if (!wbinary) 
            N.paired.rows[counter] <-
              sum(first - second > 0 & second > 0)/sum(second > 0)
          else
            N.paired.rows[counter] <-
              sum(first - second >= 0 & second > 0)/sum(second > 0)
        }
        else {
          N.paired.rows[counter] <-
            sum(first + comm[j, ] == 2)/rfill[j]
        }
      }
    }
    counter <- 0
    for (i in 1:(nc - 1)) {
      first <- comm[, i]
      for (j in (i + 1):nc) {
        counter <- counter + 1
        if (cfill[i] <= cfill[j] || any(cfill[c(i, j)] == 0)) 
          next
        if (weighted) {
          second <- comm[, j]
          if (!wbinary)
            N.paired.cols[counter] <-
              sum(first - second > 0 & second > 0)/sum(second > 0)
          else
            N.paired.cols[counter] <-
              sum(first - second >= 0 & second > 0)/sum(second > 0)
        }
        else {
          N.paired.cols[counter] <-
            sum(first + comm[, j] == 2)/cfill[j]
        }
      }
    }
    N.columns <- mean(N.paired.cols) * 100
    N.rows <- mean(N.paired.rows) * 100
    NODF <- (sum(c(N.paired.rows, N.paired.cols)) * 100)/
      ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
    out <- list(comm = comm, fill = fill,
                statistic = c(N.columns = N.columns, N.rows = N.rows, NODF = NODF))
    class(out) <- "nestednodf"
    out
  }



master = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/master4.24.csv")

master.nest = master[ c(39:70) ]
master.nest = na.omit(master.nest)
names(master.spp)

library(vegan)
worm.final.nest = oecosimu(comm = master.nest, nestfun = nestednodf, method = "swap", 
                              nsimul = 1000, alternative = c("two.sided"))

epigeic = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/Epigeic.csv")
view(epigeic)
epigeic = epigeic[ c(39:48) ]
epigeic = na.omit(epigeic)
names(master.spp)

epigeic.nest = oecosimu(comm=epigeic, nestfun = nestednodf, method = "swap", 
                        nsimul = 1000, alternative = c("two.sided"))



### Correlation anaylsis
worms.cor = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/worms.4.22.csv")
names(worms.cor)
r.year = rank(worms.cor$sample.year)
r.latitude = rank(worms.cor$latitude)

max.pack = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/max.pack.sheet.csv")
rm.years = rank(max.pack$sample.year)

rm.latitude = rank(max.pack$latitude)

max.spp = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/max.spp.csv")



cor.test(max.pack$max.pack.site, rm.years, alternative = "two.sided", 
         method = "spearman", conf.level = 0.95)
cor.test(max.pack$max.pack.site, rm.latitude, alternative = "two.sided", 
         method = "spearman", conf.level = 0.95)
names(max.pack)
r.func = rank(max.spp$func.grp)
r.func

cor.test(max.spp$max.pack.spp, r.func, alternative = "two.sided", 
         method = "spearman", conf.level = 0.95)




##visualization of max pack matrix
master = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/master4.24.csv")

master.spp = master[ c(39:68) ]
master.spp = na.omit(master.spp)
names(master.spp)

#master.spp = as.numeric(master.spp, header = TRUE)
# Not working
w = vegan::nestednodf(master.spp, order = TRUE, weighted = FALSE, wbinary = FALSE)

## S3 method for class 'nestednodf'

spp.names = c("Dendrobaena octaedra", "Lumbricus spp.", "Apporectodea spp.", 
              "Dendrobaena spp.", "Lumbricus rubellus", 
              "Lumbricus terrestris", "Octalasion tyrtaeum", "Octalasion spp.",
              "Aporrtectodea caliginosa", "Aporrectodea rosea", 
              "Bimastos rubidus", "Aporrectodea turgida", "Lumbricus castaneus",
              "Dendrodrilus spp.", "Octalasion cyaneum", "Amynthas agrestis",
              "Eisenia teteadra", "Amynthas spp.", "Endogeic spp.",
               "Allolobophora cholorotica","Amynthas tokioensis",
              "Eisenia fetida", "Eisenia hortensis", "Aporrectodea longa",
              "Eisenia lucens", "Bimastos parvus",
              "Epigeic spp.", "Epi-endogeic spp.", "Anecic spp.")
spp.names


e = vegan::nestednodf(epigeic, order = TRUE, weighted = FALSE, wbinary = FALSE)



############## FINAL CLEANED CODE FOR EARTHWORM ANALYSIS
##############################################################
master = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/master.csv")
view(master)
  
### Create a matrix for nestedness that discludes habitat and other variables, (just spp.)
# and cut out "unknown" spp. category.
master.nest = master[ c(39:67) ]
master.nest = na.omit(master.nest)
names(master.spp)
view(master.nest)

master.nest2 = master.nest[c(1:25)]
  
#Run nested NODF test for whole matrix
library(vegan)
worm.final.nest.spponly = oecosimu(comm = master.nest2, nestfun = nestednodf, method = "swap", 
                             nsimul = 1000, alternative = c("two.sided"))
  
# Visualize maximally packed matrix using vegan::nestednodf function
w = vegan::nestednodf(master.nest2, order = TRUE, weighted = FALSE, wbinary = FALSE)

#Create plot using nestednoddf "w" output

?par
par(mar = c(6.5, 5, 10, 4.1))

plot(w, col = "blue", names = TRUE, xlab = "", 
     ylab = "Site Order by Decreasing Fill", yaxt = "n", cex.lab = 2)
abline(a = 0, b = 1)


spp.names1 = as.vector(spp.names)




?axis


#Adjust x axis labels

axis(1, seq(1:30))
axis(1, at = 1:30, labels = spp.names, las = 3)
  
spp.names = (c("Dendrobaena octaedra", "Lumbricus spp.", "Apporectodea spp.", 
                "Dendrobaena spp.", "Lumbricus rubellus", 
                "Lumbricus terrestris", "Octalasion tyrtaeum", "Octalasion spp.",
                "Aporrtectodea caliginosa", "Aporrectodea rosea", 
                "Bimastos rubidus", "Aporrectodea turgida", "Lumbricus castaneus",
                "Dendrodrilus spp.", "Octalasion cyaneum", "Amynthas agrestis",
                "Eisenia teteadra", "Amynthas spp.", "Endogeic spp.",
                "Allolobophora cholorotica","Amynthas tokioensis",
                "Eisenia fetida", "Eisenia hortensis", "Aporrectodea longa",
                "Eisenia lucens", "Bimastos parvus",
                "Epigeic spp.", "Epi-endogeic spp.", "Anecic spp."))
spp.names = as.vector(spp.names)

max.pack = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/max.pack.sheet.csv")
max.pack1 = max.pack[c(44:74)]
names(max.pack1)

max.pack2 = as.matrix(max.pack1, header = FALSE)
view(max.pack2)


max.pack2 = tibble(max.pack2)

list.max = as.list(max.pack2)
head(list.max)

  
#Run Spearman correlation analysis on maximally packed sites by latitude, sample year, potentially
# temp, precip, eleveation
  
#Upload "max.pack.sheet.csv"

max.pack = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/max.pack.sheet.csv")
  
#Use max.pack.sheet$max.pack.site
  
#Rank transform environmental variables
r.years = rank(max.pack$sample.year)

r.latitude = rank(max.pack$latitude)
  
#r.temp
#r.precip
#r.elev

#Spearman Correlation
cor.test(max.pack$max.pack.site, r.years, alternative = "two.sided", 
           method = "spearman", conf.level = 0.95)
cor.test(max.pack$max.pack.site, r.latitude, alternative = "two.sided", 
           method = "spearman", conf.level = 0.95)
  
###Waiting on these data from GIS download
cor.test(max.pack$max.pack.site, r.temp, alternative = "two.sided", 
           method = "spearman", conf.level = 0.95)
cor.test(max.pack$max.pack.site, r.precip, alternative = "two.sided", 
           method = "spearman", conf.level = 0.95)
cor.test(max.pack$max.pack.site, r.elev, alternative = "two.sided", 
           method = "spearman", conf.level = 0.95)
  
#Run Spearman correlation analysis on maximally packed species by functional group

#upload "max.spp.csv" file
max.spp = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/max.spp.csv")

#Functional Group Codes, by depth into soil:
#Epigeic = 1
#Epi-endogeic = 2
#Endogeic = 3
#Anecic = 4
  
  
#Rank transform functional group list
r.func = rank(max.spp$func.grp)
  
#Spearman Correlation
cor.test(max.spp$max.pack.spp, r.func, alternative = "two.sided", 
         method = "spearman", conf.level = 0.95)
  
###Nested NODF test by functional group
  
#Upload "functional.master.csv" file (will not include genera and unknowns)
function.master = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/function.master.csv")

#Cut out all but spp. data
function.master = function.master[ c(67:70) ]
function.master[-1] = as.integer(function.master[-1] != 0)
function.master = na.omit(function.master)


#Run vegan::oecosimu
function.nest = oecosimu(comm=function.master, nestfun = nestednodf, method = "swap", 
                        nsimul = 1000, alternative = c("two.sided"))
  
##Testing for nestedness by functional group

#plot
names(function.master)

colnames(function.master) <- c("Anecic", "Endogeic", "Epigeic", "Epi-endogeic")
e = vegan::nestednodf(function.master, order = TRUE, weighted = FALSE, wbinary = FALSE)

par(mar = c(6.5, 5, 10, 4.1))


plot(w, col = "blue", names = TRUE, xlab = "", 
     ylab = "Site Order by Decreasing Fill", yaxt = "n", cex.lab = 2)
abline(a = 0, b = 1)

plot(e, col = "blue", names = TRUE, xlab = "", 
     ylab = "Site Order by Decreasing Fill", yaxt = "n", cex.lab = 2)
abline(a = 0, b = 1)
#axis(side = 3, at = 1:30, line = 1, pos = c(1,-20), labels = spp.names, 
#las = 3, hadj = 0, padj = -37)

#Fox rule evaluation

#categorized sites by meeting rule or not, by hand
# conservative method = counted epi-endogeic as separate group
# non conservative methods = did not count epi-endogeic as separte group

fox = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/Fox.total.csv")
view(fox)
fox = na.omit(fox)
t.test(fox$Yes, fox$No, var.equal = TRUE)

var.test(fox$Yes, fox$No)

fox.matrix = matrix(c(114, 88, 101, 101), nrow = 2, byrow = TRUE)
fox.matrix

row.names(fox.matrix) = c("observed", "expected")
colnames(fox.matrix) = c("Yes", "No")

fox.matrix
x = (c(114,88))
y = (c(101,101))

observed = c(114, 88)
expected = c(101, 101)

chi.test = chisq.test(fox.matrix)

summary(chi.test)
??chisq.test

##Welch two-sample t-test says Fox rule happens more often than not

###Nested NODF test for only epigeic spp.
  
#Upload "Epigeic.csv" file
epigeic = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/Epigeic.csv")
  
#Cut out all but spp. data
epigeic = read.csv("~/Desktop/Earthworm GIS/EarthwormStats/Epigeic.csv")
epigeic = epigeic[ c(39:48) ]
epigeic = na.omit(epigeic)

#Run vegan::oecosimu
epigeic.nest = oecosimu(comm=epigeic, nestfun = nestednodf, method = "swap", 
                        nsimul = 1000, alternative = c("two.sided"))

##Are epigeic species nested? Expect yes but stats will be weak.
epigeic.nest
#No, N.columns is 17.83937, p-value = 0.4727527