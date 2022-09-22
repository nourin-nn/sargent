test_that("test sargentAnnotation multiple gene sets", {
  x <- get(data("sargentDemoData"))
  srgnt <- sargentAnnotation(gex=x$gex,
                             gene.sets=x$gene.sets,
                             adjacent.mtx=x$adjacent.mtx)
  srgnt_tst <- setNames(c(19, 16, 25, 23, 17), 
                        c("B", "Endo", "Epit", "MPh", "TNK"))
  expect_equal(print(srgnt), srgnt_tst)
})


test_that("test sargentAnnotation single gene sts", {
  x <- get(data("sargentDemoData"))
  srgnt <- sargentAnnotation(gex=x$gex,
                             gene.sets=list("B"=x$gene.sets$B), 
                             score.threshold=0.25)
  srgnt_tst <- setNames(c(22, 78), 
                        c("B", "TBD"))
  expect_equal(print(srgnt), srgnt_tst)
})