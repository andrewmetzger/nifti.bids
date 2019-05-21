label.summary <- function(data.nii,
                          label.nii,
                          label.key,
                          scratch.dir = NULL,
                          save.vol = TRUE,
                          save.stats = c("voxels", "mean", "median", "sd", "mad", "skew", "kurtosis", "975", "025"),
                          overwrite.sub = TRUE,
                          save.dir=NULL,
                          save.prefix=NULL,
                          verbose.calc=FALSE) {

  # Check inputs
  stopifnot(all(file.exists(data.nii, label.nii, label.key)))

  # Retrieve input params
  temp <- unlist(strsplit(data.nii, split="/"))
  top.dir <- paste(temp[1:(which(temp=="derivatives")-1)], collapse="/")

  if (is.null(scratch.dir)) {
    scratch.dir <- paste(top.dir, "scratch", sep="/")
  }
  if (!dir.exists(scratch.dir)) { dir.create(scratch.dir) }
  # unzip niftis
  if (file_ext(data.nii) == "gz") {
    datatemp <- data.nii
    data.nii <- sprintf("%s/%s", scratch.dir, basename(file_path_sans_ext(data.nii)))
    gunzip(filename=datatemp,
           destname=data.nii, remove=FALSE, overwrite=TRUE)
  } else if (file_ext(data.nii) != "nii") {
    stop("Unknown file type for data input")
  }
  if (file_ext(label.nii) == "gz") {
    datatemp <- label.nii
    label.nii <- sprintf("%s/%s", scratch.dir, basename(file_path_sans_ext(label.nii)))
    gunzip(filename=datatemp,
           destname=label.nii, remove=FALSE, overwrite=TRUE)
  } else if (file_ext(label.nii) != "nii") {
    stop("Unknown file type for data input")
  }

  if (is.null(save.dir)) { save.dir <- paste(top.dir, "summary", sep="/") }
  if (!dir.exists(save.dir)) { dir.create(save.dir) }

  # retrieve file name and path info
  temp <- unlist(strsplit(unlist(strsplit(label.nii, split="/scratch/"))[1], split="/"))
  researcher <- paste0(temp[-length(temp)], collapse="/")
  project <- temp[length(temp)]

  subject <- unlist(strsplit(unlist(strsplit(data.nii, split="sub-"))[2], split="_"))[1]
  session <- unlist(strsplit(unlist(strsplit(data.nii, split="ses-"))[2], split="_"))[1]

  temp <- unlist(strsplit(unlist(strsplit(label.nii, split="-")), "/"))
  temp <- temp[length(temp)]
  label.name <- unlist(strsplit(temp[length(temp)], "[.]"))[1]

  if (is.null(save.prefix)) {
    save.prefix <- sprintf("%s_%s", project, label.name)
  }

  keyf <- read.delim(label.key, sep="\t", stringsAsFactors = FALSE)
  num.rows <- 0
  out.names <- character()
  if (save.vol) {
    out.names <- "volume"
  }
  if (!is.logical(save.stats) & !is.null(save.stats)) {
    out.names <- c(out.names, save.stats)
    data.type <- unlist(strsplit(unlist(strsplit(unlist(strsplit(data.nii, "/")), "[.]")), "_"))
    data.type <- data.type[length(data.type)-1]
  }
  num.rows <- length(out.names)

  pf <- data.frame(subject = rep(subject, num.rows),
                   session = rep(session, num.rows),
                   value = out.names,
                   matrix(0, nrow=num.rows, ncol=nrow(keyf)),
                   stringsAsFactors = FALSE)
  colnames(pf) <- c("subject", "session", "value", keyf$roi)

  mod.vals <- read.nii.volume(data.nii, 1)
  labels <- read.nii.volume(label.nii, 1)
  vxl.size <- prod(unlist(nii.hdr(data.nii, "pixdim"))[2:4])

  for (i in 1:nrow(keyf)) {
    roi.values <- as.numeric(unlist(strsplit(keyf$value[i], split=",")))
    roi.idx <- arrayInd(which(labels %in% roi.values), .dim=dim(labels))
    for (j in out.names) {
      pf[which(pf$value == j), i+3] <- switch(j,
        `volume` = nrow(roi.idx) * vxl.size,
        `voxels` = nrow(roi.idx),
        `mean` = mean(mod.vals[roi.idx], na.rm=TRUE),
        `median` = median(mod.vals[roi.idx], na.rm=TRUE),
        `sd` = sd(mod.vals[roi.idx], na.rm=TRUE),
        `mad` = mad(mod.vals[roi.idx], na.rm=TRUE),
        `skew` = skewness(mod.vals[roi.idx], na.rm=TRUE),
        `kurtosis` = kurtosis(mod.vals[roi.idx], na.rm=TRUE),
        `975` = quantile(mod.vals[roi.idx], as.numeric(j)/1000, na.rm=TRUE),
        `025` = quantile(mod.vals[roi.idx], as.numeric(j)/1000, na.rm=TRUE))

      if (verbose.calc) {
        print(sprintf("calculating %s %s", keyf$roi[i], j))
      }
    }
  }

  # find existing file
  if (save.vol) {
    fls <- list.files(path = save.dir,
                      pattern="^.*volumes.*.csv$",
                      full.names = TRUE)
    if (length(fls) == 0) {
      write.table(x = pf[1, ],
                  file=sprintf("%s/%s_volumes_%s.csv", save.dir, save.prefix,
                               format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                  quote=FALSE, row.names = FALSE, col.names = TRUE, sep=",")
    } else {
      if (length(fls >= 1)) { newest <- fls[length(fls)] }
      tf <- read.csv(newest, stringsAsFactors = FALSE)
      if (overwrite.sub) {
        sub.match <- which((tf$subject == pf$subject[1]) * (tf$session == pf$session[1])==1)
        if (length(sub.match) != 0) {
          tf[sub.match, ] <- pf[1, ]
        } else {
          tf <- rbind(tf,pf[1, ])
        }
      } else {
        tf <- rbind(tf,pf[1, ])
      }
      write.table(x = tf,
                  file=sprintf("%s/%s_volumes_%s.csv", save.dir, save.prefix,
                               format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                  quote=FALSE, row.names = FALSE, col.names = TRUE, sep=",")
      invisible(file.remove(newest))
    }
  }

  if (!is.logical(save.stats) & !is.null(save.stats)) {
    fls <- list.files(path = save.dir,
                      pattern=sprintf("^.*stats-%s.*.csv$", data.type),
                      full.names = TRUE)
    if (length(fls) == 0) {
      write.table(x = pf[-1, ],
                  file=sprintf("%s/%s_stats-%s_%s.csv", save.dir, save.prefix, data.type,
                               format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                  quote=FALSE, row.names = FALSE, col.names = TRUE, sep=",")
    } else {
      if (length(fls >= 1)) { newest <- fls[length(fls)] }
      tf <- read.csv(newest, stringsAsFactors = FALSE)
      if (overwrite.sub) {
        sub.match <- which((tf$subject == pf$subject[1]) * (tf$session == pf$session[1])==1)
        if (length(sub.match) != 0) {
          tf[sub.match, ] <- pf[-1, ]
        } else {
          tf <- rbind(tf,pf[-1, ])
        }
      } else {
        tf <- rbind(tf,pf[-1, ])
      }
      write.table(x = tf,
                  file=sprintf("%s/%s_stats-%s_%s.csv", save.dir, save.prefix, data.type,
                               format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                  quote=FALSE, row.names = FALSE, col.names = TRUE, sep=",")
      invisible(file.remove(fls))
    }
  }
  invisible(file.remove(data.nii))
  invisible(file.remove(label.nii))
}
