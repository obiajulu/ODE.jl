## Non-stiff fixed step solvers
nonstiff_fixedstep= [
           ODE.ode1,
           ODE.ode2_midpoint,
           ODE.ode2_heun,
           ODE.ode4,
           ODE.ode4ms,
           ODE.ode5ms
           ]

## Non-stiff fixed step solvers
nonstiff_adaptive=[
#          ODE.ode21, # this fails on Travis with 0.4?! TODO revert once fixed.
           ODE.ode23,
           ODE.ode45,
           ODE.ode45_dp,
           ODE.ode45_fe,
           ODE.ode78
           ]
# Stiff fixed-step solvers
stiff_fixedstep=[
           ODE.ode4s_s,
           ODE.ode4s_kr
           ]
#Stiff adaptive solvers
stiff_adaptive = [
           ODE.ode23s
           ]

#All of ODE.jl solvers by type
all_solvers = [nonstiff_fixedstep,
              nonstiff_adaptive,
              stiff_fixedstep,
              stiff_adaptive,
              ]

## minimal function export list
minimal_solvers = []# adaptive non-stiff:
              #export ode23, ode45, ode78
              # non-adaptive non-stiff:
              #export ode4, ode4ms
              # adaptive stiff:
              #export ode23s
              # non-adaptive stiff:
              #export ode4s
