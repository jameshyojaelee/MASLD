ex <- parse("scripts/s3a_replication.R"); t <- vapply(ex, function(e) paste(deparse(e), collapse = " "), "")
i <- grep("RB|NULLS|null_res|gc[(]|Null fits", t); for (k in i) cat(k, substr(t[k], 1, 90), "\n")
