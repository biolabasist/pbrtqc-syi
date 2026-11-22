# Builds the de-identified analysis files from the laboratory information system
# exports. Patient-level data are not distributed with this repository.
#
# Inputs (one workbook per test): adult results of the central laboratory
# (two Beckman Coulter AU5800 analyzers), 1 January to 31 December 2023, after
# exclusion of requests from emergency, intensive care, dialysis and oncology
# units; and the 2023 internal quality control export.
#
# The LIS exports (one workbook per test, named <test>_merkez_eriskin.xlsx) and
# the IQC export are read from PBRTQC_LIS and PBRTQC_IQC.
# Run from the repository root: Rscript scripts/01_prepare_data.R

source("R/setup.R")
library(openxlsx)

lis_dir <- Sys.getenv("PBRTQC_LIS", "lis_export")
iqc_file <- Sys.getenv("PBRTQC_IQC", "lis_export/iqc_2023.xlsx")

lis_names <- c(albumin = "albumin", alp = "alp", alt = "alt", ast = "ast",
               calcium = "kalsiyumtotal", chloride = "klor", hdl = "hdl", ldl = "ldl",
               cholesterol = "totalkolesterol", creatinine = "kreatinin", ggt = "ggt",
               glucose = "glukozaclik", iron = "demir", ldh = "ldh", magnesium = "magnezyum",
               phosphate = "fosforinorganik", potassium = "potasyum", sodium = "sodyum",
               bilirubin = "bilirubintotal", protein = "proteintotal",
               triglycerides = "trigliserid", urea = "ureazotu", uric_acid = "urikasit")

excluded_units <- "onkoloji|diyaliz|yogun|acil"

read_test <- function(analyte) {
  f <- file.path(lis_dir, sprintf("%s_merkez_eriskin.xlsx", lis_names[[analyte]]))
  d <- read.xlsx(f)
  unit <- tolower(paste(d$bolum, d$bolbirim))
  d <- d[!grepl(excluded_units, unit), ]
  data.frame(analyte = analyte,
             time = as.POSIXct(d$sonuctarih, tz = "UTC"),
             value = as.numeric(d$sonuc),
             analyzer = as.integer(sub(".*AU5800 ([12])$", "\\1", d$cihaz)),
             stringsAsFactors = FALSE)
}

stream <- do.call(rbind, lapply(analytes$analyte, read_test))
stream <- stream[order(stream$analyte, stream$time), ]
saveRDS(stream, file.path(data_dir, "stream_2023.rds"))

iqc <- read.xlsx(iqc_file, sheet = 1)
hdr <- unlist(iqc[1, ])
names(iqc) <- make.unique(ifelse(is.na(hdr), names(iqc), hdr))
iqc <- iqc[-1, ]
iqc <- iqc[grepl("\\.xls$", iqc$Source.Name), ]
src <- sub("\\.xls$", "", iqc$Source.Name)
test <- trimws(sub(" ?[0-9]+(m[12])?$", "", src))
test[grepl("ürikasit", src)] <- "ürikasit"
test[test == "kretinin"] <- "kreatinin"
iqc <- data.frame(iqc_name = test,
                  analyzer = ifelse(grepl("AU5800 1$", iqc$`Cihaz Adı`), 1L, 2L),
                  module = ifelse(grepl("m2", src, fixed = TRUE), 2L, 1L),
                  level = ifelse(iqc$`Kontrol Tipi` == "Düşük", 1L, 2L),
                  time = convertToDateTime(as.numeric(iqc$`Sonuç Tarihi`), tz = "UTC"),
                  value = as.numeric(iqc$Sonuç),
                  westgard = iqc$Westgard,
                  lot = iqc$`Kontrol Lot No`,
                  valid = iqc$`Geçerli Sonuç` == "Var",
                  stringsAsFactors = FALSE)
iqc <- iqc[!is.na(iqc$valid) & iqc$valid & !is.na(iqc$value), ]
saveRDS(iqc, file.path(data_dir, "iqc_2023.rds"))
