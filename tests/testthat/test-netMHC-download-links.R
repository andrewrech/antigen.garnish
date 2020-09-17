testthat::test_that("netMHC-download", {
  links <- c(
    "http://www.cbs.dtu.dk/services/NetMHC-4.0/data.tar.gz",
    "http://www.cbs.dtu.dk/services/NetMHCII-2.3/data.Linux.tar.gz",
    "http://www.cbs.dtu.dk/services/NetMHCIIpan/data.tar.gz",
    "http://www.cbs.dtu.dk/services/NetMHCpan-4.1/data.tar.gz"
  )
  for (i in links) {
    cmd <- paste0("curl --head --silent --fail ", i, "> /dev/null 2>&1")
    print(cmd)
    exit <- system(cmd)
    testthat::expect_equal(exit, 0)
  }
})
