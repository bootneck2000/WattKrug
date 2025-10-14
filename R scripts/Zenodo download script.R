library(zen4R)
library(httr)

# Your Zenodo record ID
record_id <- "5040114"

# Get Zenodo record metadata
zenodo <- ZenodoManager$new()
record <- zenodo$getRecordById(record_id)

# Create download directory if it doesn't exist
dir.create("zenodo_downloads", showWarnings = FALSE)

for (i in record$files) {
  file_name <- record$files[[1]]$filename[1]
  file_url <-  paste0("https://zenodo.org/api/records/", record_id, "/files/", file_name, "/content")
  
  # URLencode only the file name part, not the whole URL
  dest_path <- file.path("zenodo_downloads", file_name)
  cat("Downloading", file_name, "from", file_url, "...\n")
  resp <- GET(file_url, write_disk(dest_path, overwrite = TRUE), add_headers(`User-Agent` = "R"))
  stop_for_status(resp)
  cat("Downloaded:", dest_path, "\n")
}

cat("All files downloaded!\n")
#then move WD to the downloaded files.
setwd("./zenodo_downloads")
