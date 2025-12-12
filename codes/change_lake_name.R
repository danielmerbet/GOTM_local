all_lakes = c("Allequash", "Alqueva", "Annie", "Arendsee", "Argyle", "Biel", 
              "BigMuskellunge", "BlackOak", "Bosumtwi", "Bryrup", "BurleyGriffin", 
              "Chao", "Crystal", "CrystalBog", "Delavan", "Dickie", "Eagle", 
              "Ekoln", "Erken", "EsthwaiteWater", "FallingCreek", "Feeagh", 
              "Fish", "GreatPond", "Green", "Harp", "Hassel", "Hulun", "Kilpisjarvi", 
              "Kinneret", "Kivu", "Klicava", "Kuivajarvi", "Langtjern", "Laramie", 
              "LowerLakeZurich", "Mendota", "Monona", "Mozhaysk", "MtBold", 
              "Muggelsee", "Murten", "Neuchatel", "Ngoring", "NohipaloMustjarv", 
              "NohipaloValgejarv", "Okauchee", "Paajarvi", "Rappbode", "Rappbodep", 
              "Rimov", "Rotorua", "Sammamish", "Sau", "Scharmutzel", "Sparkling", 
              "Stechlin", "Sunapee", "Tahoe", "Taihu", "Tarawera", "Thun", 
              "Toolik", "Trout", "TroutBog", "TwoSisters", "Vendyurskoe", "Vortsjarv", 
              "Washington", "Windermere", "Wingra", "Zlutice", "Zurich")
tolower(all_lakes)
cal <- "calibrated"
for (e in c("historical", "ssp126", "ssp370", "ssp585", "picontrol")){
  for (g in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
    a <- read.csv(paste0("~/Documents/isi_cal/", cal,"/annual/",e,
                         "/GOTM_",e,"_", tolower(g),"_", cal, ".csv"))
    
    a$lake <- sapply(a$lake, function(x) {
      match_name <- all_lakes[match(tolower(x), tolower(all_lakes))]
      if (is.na(match_name)) warning(paste("Lakename not found:", x))
      match_name
    })
    
    write.csv(a, paste0("~/Documents/isi_cal/",cal,"/annual_ok/",e,
                    "/GOTM_",e,"_", tolower(g),"_", cal, ".csv"), 
              quote = F, row.names = F)
  }
}
