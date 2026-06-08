## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(sparsediff)
engine_version()

## ----build--------------------------------------------------------------------
n   <- 3L
x   <- sd_variable(d1 = n, d2 = 1L, var_id = 0L, n_vars = n)  # 0-based var offset
obj <- sd_sum(sd_exp(x), axis = -1L)                          # sum(exp(x))  (axis -1 = all)
g1  <- sd_sum(x, axis = -1L)                                  # sum(x)

## ----problem------------------------------------------------------------------
prob <- sd_problem(objective = obj, constraints = list(g1), verbose = FALSE)

# initialise the (structural) derivative layout once
sd_init_derivatives(prob)
sd_init_jacobian_coo(prob)
sd_init_hessian_coo(prob)

## ----evaluate-----------------------------------------------------------------
u <- c(0, 0.5, 1)

sd_objective_forward(prob, u)   # value of sum(exp(x))
sd_gradient(prob)               # gradient = exp(u)

sd_constraint_forward(prob, u)  # value of sum(x)

## ----jacobian-----------------------------------------------------------------
sd_jacobian_sparsity(prob)      # rows / cols / nrow / ncol
sd_jacobian_values(prob)        # d/dx sum(x) = (1, 1, 1)

## ----hessian------------------------------------------------------------------
sd_hessian_sparsity(prob)
sd_hessian_values(prob, obj_w = 1, w = 0)

## ----parameters---------------------------------------------------------------
p    <- sd_parameter(d1 = 1L, d2 = 1L, param_id = 0L, n_vars = n, values = 2)
obj2 <- sd_sum(sd_scalar_mult(p, sd_exp(x)), axis = -1L)      # theta * sum(exp(x))
prob2 <- sd_problem(obj2, constraints = list(), verbose = FALSE)
sd_register_params(prob2, list(p))
sd_init_derivatives(prob2)

sd_objective_forward(prob2, u)        # theta = 2
sd_update_params(prob2, 3)            # theta -> 3
sd_objective_forward(prob2, u)        # re-evaluated, no rebuild

