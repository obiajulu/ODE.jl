# ODE23S  Solve stiff systems based on a modified Rosenbrock triple
# (also used by MATLAB's ODE23s); see Sec. 4.1 in
#
# [SR97] L.F. Shampine and M.W. Reichelt: "The MATLAB ODE Suite," SIAM Journal on Scientific Computing, Vol. 18, 1997, pp. 1–22
#
# supports keywords: points = :all | :specified (using dense output)
#                    jacobian = G(t,y)::Function | nothing (FD)

# Internal

immutable ModifiedRosenbrockStepper{T<:Number} <: AbstractStepper
    d :: T
    e32 :: T

    function ModifiedRosenbrockStepper()
        d   = T(1/(2 + sqrt(2)))
        e32 = T(6 + sqrt(2))
        new(d,e32)
    end
end


# TODO: is this correct?
order(ModifiedRosenbrockStepper) = 2


# define the set of ODE problems with which this stepper can work
solve(ode :: ExplicitODE, stepper :: ModifiedRosenbrockStepper, options :: Options) =
    Solver{ModifiedRosenbrockStepper}(ode,stepper,options)


# lower level interface (iterator)

"""
The state for the Rosenbrock stepper

- step:  Last successful step
- F1,F2: Work arrays for storing the intermediate values of y'
- J:     Jacobian
- iters: Number of successful steps made

"""
type RosenbrockState{T,S} <: AbstractState
    step ::Step{T,S}
    dt   ::T
    F1   ::S
    F2   ::S
    J # :: ?
    iters::Int
end


# for debugging
function show(io::IO, state :: RosenbrockState)
    show(io,state.step)
    println("dt =$(state.dt)")
    println("F1 =$(state.F1)")
    println("F2 =$(state.F2)")
    println("J  =$(state.J)")
end


function start(s :: Solver{ModifiedRosenbrockStepper})
    t  = s.ode.t0
    dt = s.options.initstep
    y  = s.ode.y0
    dy = zero(y)

    J  = Array(eltype(y),length(y),length(y))

    step  = Step(t,deepcopy(y),deepcopy(dy))
    state = RosenbrockState(step,
                            dt,
                            zero(y), # F1
                            zero(y), # F2
                            J,       # J
                            0)       # iters
    # initialize the derivative and the Jacobian
    s.ode.F!(t,y,step.dy)
    s.ode.jac!(t,y,state.J)

    return state
end


function next(s     :: Solver{ModifiedRosenbrockStepper},
              state :: RosenbrockState)

    stepper = s.stepper
    ode     = s.ode
    step    = state.step
    opts    = s.options

    F1, F2, J = state.F1, state.F2, state.J

    t, dt, y, dy = step.t, state.dt, step.y, step.dy
    # F!, jac! = ode.F!, ode.jac!
    d, e32 = stepper.d, stepper.e32

    F0 = dy

    while true

        state.iters += 1
        if state.iters > s.options.maxiters
            return ((step.t,step.y), state)
        end

        # trim the step size to match the bounds of integration
        dt = min(s.options.tstop-t,dt)

        # TODO: this should go to a specialized function for type stabilty sake
        # maybe make W a part of ExplicitODE?  Same for tder below?
        if size(J,1) == 1
            W = one(J) - dt*d*J
        else
            # note: if there is a mass matrix M on the lhs of the ODE, i.e.,
            #   M * dy/dt = F(t,y)
            # we can simply replace eye(J) by M in the following expression
            # (see Sec. 5 in [SR97])

            W = lufact( eye(J) - dt*d*J )
        end

        # Approximate time-derivative of F, we are using F1 as a
        # temporary array
        ode.F!(t+dt/100,y,F1)
        tder = 100*d*(F1-F0)

        # modified Rosenbrock formula
        # TODO: allocate some temporary space for these variables
        k1 = W \ (F0 + tder)
        ode.F!(t+dt/2, y+dt*k1/2, F1)
        k2 = W \ (F1 - k1) + k1
        ynew = y + dt*k2
        ode.F!(t+dt,   ynew,      F2)
        k3 = W \ (F2 - e32*(k2 - F1) - 2*(k1 - F0) + tder )

        delta = max(opts.reltol*max(opts.norm(y),
                                    opts.norm(ynew)),
                    opts.abstol) # allowable error

        err = (dt/6)*opts.norm(k1 - 2*k2 + k3)/delta # error estimate

        # upon a failed step decrease the step size
        dtnew = min(opts.maxstep,
                    dt*0.8*err^(-1/3) )

        # check if the new solution is acceptable
        if err <= 1

            # update the state and return
            step.t     = t+dt
            state.dt   = dtnew
            step.y[:]  = ynew
            step.dy[:] = F2
            ode.jac!(step.t,step.y,J)

            return ((step.t,step.y), state)
        else
            # continue with the decreased time step
            dt = dtnew
        end

    end

    return tout, yout

end
