#=
    (1) done()
    (2) trialstep!()
    (3) errorcontol!()
    (4) odercontrol!()
    (5) rollback!()
    (6) status()
=#

<<<<<<< 05b0e332ba8054e68019fd93b1815b48e5485d87
type RadauState{T,Y}
    h::T     # (proposed) next time step

    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y # derivative at t-1

    step::Int # current step number
    finished::Bool # true if last step was taken

    # work arrays
end

function radau(f, y0, tspan, order ::Integer = 5)
=======
function radau(f, y0, tspan, order::Integer = 5)
    #= setup
        state.t, state.f(t)
    =#
>>>>>>> add basic errorcontrol! function

    # Set up
    T = eltype(tspan)
    Y = typeof(y0)
    EY = eltype(Y)
    N = length(tsteps)
    dof = length(y0)
    h = hinit
    t = ode.tspan[1]
    y = deepcopy(y0)
    dy = f(t, y)
    ## previous data set to null to begin
    tpre = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)
    step = 1
    finished = false

    # Calculate c_i
    c
    # Calculate b_i
    b
    # Calculate a_ij
    a
    
    ## intialize state
    st =  RadauState{T,Y}(h,
                   t, y, dy,
                   tpre, ypre, dypre,
                   step, finished)

    # Time stepping loop
    while !done()
        stats = trialstep!(st)
        err, stats, st = errorcontrol!(st) # (1) adaptive stepsize (2) error
        if err < 1
            stats, st = ordercontrol!()
            accept = true
        else
            rollback!()
        end

        return = status()
    end
end

function done(st)
    @unpack st: h, t, tfinal
    if h < minstep || t = tfinal
        return true
    else
        return false
    end
end

<<<<<<< 05b0e332ba8054e68019fd93b1815b48e5485d87
function trialstep!(st)
    @unpack st: h, t, tfinal
    # Calculate simplified Jacobian
    #     _                                             _
    #    |M - h*a[11]*h*f(tn,yn) ... -h*a[1s]*f(tn,yn)   |
    # J= |         ⋮              ⋱           ⋮           |
    #    | - h*a[1s]*h*f(tn,yn) ...   M-h*a[ss]*f(tn,yn) |
    #    |_                                             _|
    #

    # Use Netwon interation
    #
    #   ⃗k^(i+1) = ⃗k ^ (i) - J^{-1}*G(⃗k)
    # where
    #           /       \
    #           |k_1^(i)|
    #   ⃗k^(i) = |   ⋮   |
    #           |k_s^(i)|
    #           \       /
    #
    #   and
    #
    #           /            \
    #          |f_1(k_1,…,k_s)|
    #  G(k) =  |      ⋮       |  = 0
    #          |f_s(k_1,…,k_s)|
    #           \            /
    #
    #
    #   and
    #
    #  f_i = k_i - f(t_n + c_ih, y_n + h*a_11*k1 + ⋯ + h*a_1s*k_s)
    #
    #

    # Once Netwon method converges after some m steps, estimated next step size
    #
    #   y = ypre + h ∑b_i*k_i^{m}
    #

end
=======
function errorcontrol!(st)
    @unpack st:M, h, A, J, tpre, ypre, b̂, b, c, g, order_number, fac
    γ0 = filter(λ -> imag(λ)==0, eig(inv(A)))

    for i = 1:order_number
        sum_part += (b̂[i] - b[i]) * h *f(xpre + c[i] * h, g[i])
    end

    yhat_y = γ0 * h * f(tpre, ypre) + sum_part
    err = norm( inv(M - h * λ * J) * (yhat_y) )

    hnew = fac * h * err_norm^(-1/4)
>>>>>>> add basic errorcontrol! function
