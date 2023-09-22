df1 <- data.frame(V1 = c('chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr3', 'chr3', 'chr4'),
                  V2 = c(10, 20, 98, 10, 35, 50, 60, 5),
                  V3 = c(25, 100, 101, 15, 46, 55, 90, 100))

df2 <- data.frame(V1 = c('chr1', 'chr1', 'chr2', 'chr2', 'chr2', 'chr3', 'chr4'),
                  V2 = c(95, 200, 45, 49, 55, 50, 101),
                  V3 = c(105, 205, 50, 51, 90, 100, 110))


require(data.table)
setDT(df1)
setDT(df2)

setkey(df1, V1,V2,V3)
setkey(df2, V1,V2,V3)

any_overlaps_dt = function(df1, df2) {
  olaps = foverlaps(df1, df2, mult="first", type="any", which=TRUE)
  as.integer(!is.na(olaps))
}

olaps_12 = any_overlaps_dt(df1, df2)
# [1] 0 1 1 0 1 1 1 0

olaps_21 = any_overlaps_dt(df2, df1)
# [1] 1 0 1 0 0 1 0