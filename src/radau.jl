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
using Polynomials
using ForwardDiff
using ODE
using Parameters


immutable TableauRKImplicit{Name, S, T} <: ODE.Tableau{Name, S, T}
    order::Integer # the order of the method
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauRKImplicit(order,a,b,c)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert S==length(c)==length(b)
        @assert length(b)==(order+1)/2
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
                                  [5//12  -1//12
                                   3//4    1//4],
                                  [3//4, 1//4]',
                                  [1//3, 1])

const bt_radau5 = TableauRKImplicit(:radau5,5, Rational{Int64},
                                [11/45-7*√(6)/360       37/225-169*√(6)/1800    -2/225 + √(6)/75
                                 37/225+169*√(6)/1800   11/45+7*√(6)/360        -2/225 - √(6)/75
                                 4/9-√(6)/36            4/9+√(6)/36             1/9             ],
                                [4/9-√(6)/36            4/9+√(6)/36             1/9             ]',
                                [2/5-√(6)/10,           2/5+√(6)/10,            1               ])


###########################################
# State for Radau Solver
###########################################
type RadauState{T,Y}
    tfinal ::T
    tdir :: Integer
    minstep ::Float64

    h::T     # (proposed) next time step

    t::T      # current time
    y::Y      # current solution
    dy::Y     # current derivative

    tpre::T  # time t-1
    ypre::Y  # solution at t-1
    dypre::Y # derivative at t-1

    step::Int # current step number
    finished::Bool # true if last step was taken

    btab #<:TableauRKImplicit # tableau according to stage number
    stageNum ::Integer
    M ::Array{Float64,2}

    # work arrays
end

###########################################
# Radau Solver
###########################################
function ode_radau(f, y0, tspan, order ::Integer = 5;   M = eye(length(y0),length(y0)),
                                                        minstep = abs(tspan[length(tspan)]-tspan[1])/10^12
                                                        )
    # Set up
    stageNum = Integer((order+1)/2)
    T = eltype(tspan)
    Y = typeof(y0)
    EY = eltype(Y)
    dof = length(y0)

    # TODO: h = hinit()
    tfinal = tspan[length(tspan)]
    tdir = sign(tfinal - tspan[1])
    h = .01
    t = tspan[1]
    y = deepcopy(y0)
    dy = f(t, y)

    ## previous data set to null to begin
    tpre = NaN
    ypre = zeros(y0)
    dypre = zeros(y0)

    step = 1
    finished = false

    stageNum = stageNum
    # get right tableau for stage number
    if stageNum ==3
        btab = bt_radau3
    elseif stageNum ==5
        btab = bt_radau5
    else
        @show typeof(stageNum)
        @show constRadauTableau(stageNum)
        btab = constRadauTableau(stageNum)
    end

    @show typeof(tfinal)
    ## intialize state
    st =  RadauState{T,Y}(tfinal,tdir, minstep,h,
                   t, y, dy,
                   tpre, ypre, dypre,
                   step, finished,
                   btab, stageNum, M)

    println("Good so far!")
    # Time stepping loop
    while !done_radau(st)
        stats = trialstep!(st)
        err, stats, st = errorcontrol!(st) # (1) adaptive stepsize (2) error
        if err < 1
            stats, st = ordercontrol!()
            accept = true
        else
            rollback!()
        end

        return status()
    end
end


 function f(t, y)
     # Extract the components of the y vector
     (x, v) = y

     # Our system of differential equations
     x_prime = v
     v_prime = -x

     # Return the derivatives as a vector
     [x_prime; v_prime]
 end


 # Initial condtions -- x(0) and x'(0)
 y = [0.0; 0.1]

 # Time vector going from 0 to 2*PI in 0.01 increments
 t = 0:0.1:4*pi;

ode_radau(f,y,t,7)
###########################################
# iterator like helper functions
###########################################

function done_radau(st)
    @unpack st: h, t, tfinal, minstep, tdir
    if h < minstep || tdir*t >= tdir*tfinal
        return true
    else
        return false
    end
end

"trial step for ODE with mass matrix"
function trialstep!(st)
    @unpack st: h, t,y, tfinal, btab

    # Form ⃗z from y

    # Calculate  Jacobian
    g(z) = f(t,z)
    J = ForwardDiff.jacobian(g, y)

    # Use to form simplied Netwon iteration matrix
    #     _                                             _
    #    |M - h*a[11]*h*f(tn,yn) ... -h*a[1s]*f(tn,yn)   |
    # G= |         ⋮              ⋱           ⋮           |
    #    | - h*a[1s]*h*f(tn,yn) ...   M-h*a[ss]*f(tn,yn) |
    #    |_                                             _|
    #
    I_N = eye(dof,dof)
    I_s = eye(stageNum,stageNum)
    M = rand(dof,dof)
    AoplusJ = kron(btab.a,J)
    IoplusM = kron(I_s,M)
    G =  IoplusM-h*AoplusJ

    # Use Netwon interation (TODO: use the transformation T^(-1)A^(-1)T = Λ,
    # W^k = (T^-1 ⊕ I)Z^k version of iteration)

    ## initial variables iteration
    ## Initialize
    z = zeros(dof)     #TODO: use better initial values for z
    zpre = zeros(dof)
    Δzpre = zeros(dof)
    Δz = zeros(dof)
    κ

    ## Matrices used for one round of iteration
    Ginv = inv(G)
    Ginv_block = Array{Float64,2}[Ginv[i*stageNum + [1:dof], j*stageNum+[1:dof]] for i = 0:stageNum-1, j= 0:stageNum-1]
    AoplusI_block = Array{Float64,2}[btab.a[i,j]*I_N for i=1:stageNum, j=1:stageNum]

    iterate = true
    count = 0
    while iterate
        Δz = Ginv_block*(-zpre + h*AoplusI_block*F(f,z,y,t,c,h))
        z = zpre + Δz  #   ⃗z^(k+1) = ⃗z ^ (k) - Δ⃗z^(k)

        # Stop condition for the Netwon iteration
        if (count >=1)
            Θ = norm(Δz)/norm(Δzpre)
            η = Θ/(1-Θ)
            if η*norm(Δz) <= κ*min(reltol,abstol)
                iterate = false
                break
            end
        end

        zpre = z
        Δzpre = Δz

        count+=1
    end

    # Once Netwon method converges after some m steps, estimated next step size
    #
    #   y = ypre + h ∑b_i*k_i^{m}
    #
    d = inv(btab.a)*btab.b
    ynext = ypre
    for i = 1 : stageNum
        ynext += z[i]*d[i]
    end
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

end
###########################################
# Other help functions
###########################################
function constRadauTableau(stageNum)
    # Calculate c_i, which are the zeros of
    #              s-1
    #            d       s-1        s
    #           ---    (x    * (x-1)  )
    #              s-1
    #           dx
    roots = zeros(Float64, stageNum - 1)
    append!(roots, [1 for i= 1:stageNum])
    poly = Polynomials.poly([roots;])
    for i = 1 : stageNum-1
        poly = Polynomials.polyder(poly)
    end
    C = Polynomials.roots(poly)

    ################# Calculate b_i #################
    #    s
    #   ___
    #   \                   1
    #   /   bᵢcᵢ^(m - 1) = ---      m = 1, ..., s
    #   ---                 m  ,
    #   i=1
    #
    # Construct a matrix C_meta to calculate B
    C_meta = Array(Float64, stageNum, stageNum)
    for i = 1:stageNum
        C_meta[i, :] = C .^ (i - 1)
    end

    # Construct a matrix 1 / stageNum
    B_meta = Array(Float64, stageNum, 1)
    for i = 1:stageNum
        B_meta[i, 1] = 1 / i
    end

    # Calculate b_i
    C_big = inv( C_meta )
    B = C_big * B_meta
    ################# Calculate a_ij ################
    #    s
    #   ___
    #   \                    cᵢ^m
    #   /   aₐⱼcᵢ^(m - 1) = -------     m = 1, ..., s
    #   ---                    m    ,   j = 1, ..., s
    #   j=1
    #
    # Construct matrix A
    A = Array(Float64, stageNum, stageNum)

    # Construct matrix A_meta
    A_meta = Array(Float64, stageNum, stageNum)
    for i = 1:stageNum
        for j = 1:stageNum
            A_meta[i,j] = B_meta[i] * C[j]^i
        end
    end

    # Calculate a_ij
    for i = 1:stageNum
        A[i,:] = C_big * A_meta[:,i]
    end
    order = 2*stageNum-1
    return TableauRKImplicit(symbol("radau$(order)"),order, A, B, C)
end


" Calculates the array of derivative values between t and tnext"
function F(f,z,y,t,c,h)
    return Array{Float64,1}[f(t+c[i]*h, y+z[i]) for i=1:stageNum]
end
