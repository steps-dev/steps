context('transition_matrix-class')

test_that('transition_matrix classes work', {
  
  # the types of transition_matrix
  mat <- matrix(c(.53,0,.42,0.1,0.77,0,0,0.12,0.9),nrow = 3,ncol = 3,byrow = TRUE)
  tmat <- as.transition_matrix(mat)
  
  # check as.transition_matrix won't handle a silly function
  expect_error(as.transition_matrix(function() x,
                           param = list(p = 0.5)))
  expect_error(as.transition_matrix(function(x) x,
                           param = list(p = 0.5)))
  expect_error(as.transition_matrix(function(x, y) x,
                           param = list(p = 0.5)))
  
  # compound transition_matrixs
  compound <- prob * rate
  compound_disp <- prob * disp
  compound_user <- prob * dd
  
  # check they have the right class
  expect_s3_class(prob, 'transition_matrix')
  expect_s3_class(rate, 'transition_matrix')
  expect_s3_class(disp, 'transition_matrix')
  expect_s3_class(dd, 'transition_matrix')
  expect_s3_class(compound, 'transition_matrix')
  expect_s3_class(compound_disp, 'transition_matrix')
  expect_s3_class(compound_user, 'transition_matrix')
  
  # check is.transition_matrix works on transition_matrixs
  expect_true(is.transition_matrix(prob))
  expect_true(is.transition_matrix(rate))
  expect_true(is.transition_matrix(disp))
  expect_true(is.transition_matrix(dd))
  expect_true(is.transition_matrix(compound))
  expect_true(is.transition_matrix(compound_disp))
  expect_true(is.transition_matrix(compound_user))
  
  # check is.transition_matrix works on non-transition_matrixs
  expect_false(is.transition_matrix(list()))
  expect_false(is.transition_matrix(NA))
  expect_false(is.transition_matrix(NULL))
  
  # check print.transfun works on boring transfuns
  expect_equal(capture.output(print(prob)),
               'probability transfun with expectation 0.5')
  expect_equal(capture.output(print(rate)),
               'rate transfun with expectation 3')
  expect_equal(capture.output(print(disp)),
               'dispersal transfun with expectation 1')
  expect_equal(capture.output(print(dd)),
               'user-specified probability transfun')
  expect_equal(capture.output(print(compound)),
               'compound transfun with expectation 1.5')
  expect_equal(capture.output(print(compound_disp)),
               'compound transfun with expectation 1')
  expect_equal(capture.output(print(compound_user)),
               'user-specified compound transfun')
  
  # screw with some transfuns and expect an error
  bad_prob2 <- bad_prob <- prob
  class(bad_prob) <- c('flooflah', 'transfun', 'function')
  class(bad_prob2) <- c('probability', 'rate', 'transfun', 'function')
  
  # they're still transfuns, but the internal checks should error
  expect_true(is.transfun(bad_prob))
  expect_true(is.transfun(bad_prob2))
  expect_error(pop:::transfunType(bad_prob))
  expect_error(pop:::transfunType(bad_prob2))
  
  # check that dispersals and rates are combined in the right way
  landscape <- as.landscape(list(coordinates = data.frame(x = runif(5),
                                                          y = runif(5)),
                                 area = data.frame(area = 1),
                                 population = data.frame(bees = 1),
                                 features = data.frame()[1, ]))
  
  disp <- d(3)
  p_disp1 <- p(0.5) * disp
  p_disp2 <- disp * p(0.2)
  
  # evaluate
  disp_mat <- disp(landscape)
  p_disp1_mat <- p_disp1(landscape)
  p_disp2_mat <- p_disp2(landscape)
  
  # check rowSums are (nearly) all 1, or the rate
  eps <- sqrt(.Machine$double.eps)
  expect_true(all((abs(rowSums(disp_mat) - 1)) < eps))
  expect_true(all((abs(rowSums(p_disp1_mat) - 1)) < eps))
  expect_true(all((abs(rowSums(p_disp2_mat) - 1)) < eps))
  
  # check that it errors when doing illegal things with dispersal transfuns
  expect_error(r(3) * disp)
  expect_error(tr(bee ~ bull, disp))
  
})