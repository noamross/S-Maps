# Font table!

fonttable <- read.table(header=TRUE, sep=",", stringsAsFactors=FALSE,
                        con <- textConnection('
Short,Canonical
mono,Courier
sans,Helvetica
serif,Times
,AvantGarde
,Bookman
,Helvetica-Narrow
,NewCenturySchoolbook
,Palatino
,URWGothic
,URWBookman
,NimbusMon
URWHelvetica, NimbusSan
,NimbusSanCond
,CenturySch
,URWPalladio
URWTimes,NimbusRom
'))
close(con)

fonttable$pos <- 1:nrow(fonttable)

library(reshape2)
fonttable <- melt(fonttable, id.vars="pos", measure.vars=c("Short","Canonical"),
                  variable.name="NameType", value.name="Font")

# Make a table of faces. Make sure factors are ordered correctly
facetable <- data.frame(Face = factor(c("plain","bold","italic","bold.italic"),
                                      levels = c("plain","bold","italic","bold.italic")))

fullfonts <- merge(fonttable, facetable)

library(ggplot2)
pf <- ggplot(fullfonts, aes(x=NameType, y=pos)) + 
  geom_text(aes(label=Font, family=Font, fontface=Face)) +
  facet_wrap(~ Face, ncol=2)

pf