
library(GEOquery)
library(RCurl)
library(GEOmetadb)


### Pull all 450K studies under GPL13534
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListFields(con, "gsm")
x<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534'")

blood_meta<-meta[unique(c(grep("blood|Blood",meta$source_name_ch1),
                          grep("blood|Blood",meta$title), 
                          grep("blood|Blood",meta$characteristics_ch1))) ,]
print(paste("Blood Samples: ",nrow(blood_meta), "   Blood Studies: ", length(unique(blood_meta$series_id)), sep=""))

#### ON KANDY