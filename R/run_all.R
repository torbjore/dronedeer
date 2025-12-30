files <- dir("R")
files <- paste0("R/", files[substring(files, 1, 3) == "fit"]) |> as.list()
lapply(files, function(i) try(source(i)))
