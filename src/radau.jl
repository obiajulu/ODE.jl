#=
    (1) done()
    (2) trialstep!()
    (3) errorcontol!()
    (4) odercontrol!()
    (5) rollback!()
    (6) status()
=#

###########################################
# Tableaus for implicit Runge-Kutta methods
###########################################

immutable TableauRKImplicit{Name, S, T} <: Tableau{Name, S, T}
    order::Integer # the order of the method
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauRKImplicit(order,a,b,c)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert c[1]==0
        @assert istril(a)
        @assert S==length(c)==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        @assert norm(sum(a,2)-c'',Inf)<1e-10 # consistency.
        new(order,a,b,c)
    end
end

function TableauRKImplicit{T}(name::Symbol, order::Integer,
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauRKImplicit{name,length(c),T}(order, a, b, c)
end
function TableauRKImplicit(name::Symbol, order::Integer, T::Type,
                   a::Matrix, b::Matrix, c::Vector)
    TableauRKImplicit{name,length(c),T}(order, convert(Matrix{T},a),
                                        convert(Matrix{T},b), convert(Vector{T},c) )
end

## Tableaus for implicit RK methods
const bt_radau3 = TableauRKImplicit(:radau3,3, Rational{Int64},
                                  [0  0
                                   1  0],
                                  [1//2, 1//2]',
                                  [0, 1])

const bt_radau5 = TableauRKImplicit(:radau5,5, Rational{Int64},
                                [0  0
                                 1  0],
                                [1//2, 1//2]',
                                [0, 1])

const bt_radau9 = TableauRKImplicit(:radau9,9, Rational{Int64},
                                [0  0
                                 1  0],
                                [1//2, 1//2]',
                                [0, 1])


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

function radau(f, y0, tspan, stageNum ::Integer = 5)
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
    stageNum = stageNum
    finished = false


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

function raduaTable(stageNum)
    # Calculate c_i, which are the zeros of
    #              s-1
    #            d       s-1        s
    #           ---    (x    * (x-1)  )
    #              s-1
    #           dx
    roots = zeros(stageNum-1)
    append!(roots, [1 for i= 1:stageNum])
    poly = Polynomials.poly([roots])
    for i = 1 : stageNum-1
        poly = Polynomials.polyder(poly)
    end
    @show Polynomials.roots(poly)

    # Calculate b_i
    b
    # Calculate a_ij
    a
    return TableauRKImplicit(order,a,b,c)
end

function done(st)
    @unpack st: h, t, tfinal
    if h < minstep || t = tfinal
        return true
    else
        return false
    end
end

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

    # Stop condition for the Netwon iteration

    # Once Netwon method converges after some m steps, estimated next step size
    #
    #   y = ypre + h ∑b_i*k_i^{m}
    #

end

function errorcontrol!(st)
    @unpack st:M, h, A, J, tpre, ypre, b̂, b, c, g, order_number, fac
    γ0 = filter(λ -> imag(λ)==0, eig(inv(A)))

    for i = 1:order_number
        sum_part += (b̂[i] - b[i]) * h *f(xpre + c[i] * h, g[i])
    end

    yhat_y = γ0 * h * f(tpre, ypre) + sum_part
    err = norm( inv(M - h * λ * J) * (yhat_y) )

    hnew = fac * h * err_norm^(-1/4)

    # Update state
    st.hnew = hnew

    return err, Nothing, st
end

function ordercontrol!(st)
    @unpack st:W, step, order
    
    if step < 10
        # update state
        st.order = 5
        
    else
        θk = norm(W[  ]) / norm(W[  ])
        θk_1 = norm(W[  ]) / norm(W[  ])

        Φ_k = √(θk * θk_1)
    end
