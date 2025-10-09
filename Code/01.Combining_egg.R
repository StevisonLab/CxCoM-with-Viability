
folder <- "Raw Data/Raw_Egg_Counts/"

# find files

files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
file_melika <- files[grepl("melika", tolower(basename(files)))]
file_eva    <- files[grepl("eva", tolower(basename(files)))]
# read (adjust header/change sep if needed)

melika <- read.csv(file_melika[1], stringsAsFactors = FALSE)
eva <- read.csv(file_eva[1],    stringsAsFactors = FALSE)

# ensure column names: find ID and count columns (simple)

id_col_m <- names(melika)[1]      # assume ID is first column
count_col_m <- names(melika)[2]   # assume count is second column
id_col_e <- names(eva)[1]
count_col_e <- names(eva)[2]

# rename for clarity

names(melika)[names(melika) == id_col_m]    <- "ID"
names(melika)[names(melika) == count_col_m] <- "CountMelika"
names(eva)[names(eva) == id_col_e]    <- "ID"
names(eva)[names(eva) == count_col_e] <- "CountEva"

# merge by ID (keep all IDs)

combined <- merge(melika, eva, by = "ID", all = TRUE)

# Use regular expressions to extract the parts
combined$Brood <- sub(".*vial14([A-Z])[0-9]+_validated", "\\1", combined$ID)
combined$VialID <- sub(".*vial14[A-Z]0*([0-9]+)_validated", "\\1", combined$ID)

# Remove the ID column
combined$ID <- NULL

# Optional: reorder columns (Brood, VialID first)
combined <- combined[, c( "VialID","Brood", "CountMelika", "CountEva")]

# view head and save
write.csv(combined, "Raw Data/Raw_Egg_Counts/Combined_eggcount.csv", row.names = FALSE)

