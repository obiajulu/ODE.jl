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
                                  [5//12  -1//12
                                   3//4    1//4],
                                  [3//4, 1//4]',
                                  [1//3, 1])

const bt_radau5 = TableauRKImplicit(:radau5,5, Rational{Int64},
                                [11/45 - 4*sqrt(6)/360  37/225 - 169*sqrt(6)/1800  -2/225 + sqrt(6)/75
                                11/45 - 4*sqrt(6)/360  37/225 - 169*sqrt(6)/1800  -2/225 + sqrt(6)/75
                                11/45 - 4*sqrt(6)/360  37/225 - 169*sqrt(6)/1800  -2/225 + sqrt(6)/75]',
                                [2//5- sqrt(6)/10, 1])

const bt_radau9 = TableauRKImplicit(:radau9,9, Rational{Int64},
                                [0  0
                                 1  0],
                                [1//2, 1//2]',
                                [0, 1])


###########################################
# State for Radau Solver
###########################################
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

    btab::TableauRKImplicit # tableau according to stage number

    # work arrays
end

###########################################
# Radau Solver
###########################################
function ode_radau(f, y0, tspan, stageNum ::Integer = 5)
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
    # get right tableau for stage number
    if stageNum ==3
        btab = bt_radau3
    elseif stageNum ==5
        btab = bt_radau5
    elseif stageNum == 9
        btab = bt_radau9
    else
        btab = constRadauTableau(stageNum)
    end

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
                   step, finished,
                   btab)

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

###########################################
# iterator like helper functions
###########################################

function done(st)
    @unpack st: h, t, tfinal
    if h < minstep || t = tfinal
        return true
    else
        return false
    end
end

"trial step for ODE with mass matrix"
function trialstep!(st)
    @unpack st: h, t,y, tfinal


    # Calculate simplified Jacobian if My' = f(t,y)
    #     _                                             _
    #    |M - h*a[11]*h*f(tn,yn) ... -h*a[1s]*f(tn,yn)   |
    # G= |         ⋮              ⋱           ⋮           |
    #    | - h*a[1s]*h*f(tn,yn) ...   M-h*a[ss]*f(tn,yn) |
    #    |_                                             _|
    #
    g(z) = f(t,z)
    J = ForwardDiff.jacobian(g, y)
    I = eye(stageNum,stageNum)

    #AoplusJ = [btab.a[i,j]*J for i=1:stageNum, j=1:stageNum]
    AoplusJ=zeros(stageNum*dof,stageNum*dof)
    for i=1:stageNum
        for j = 1:stageNum
            for l = 1:dof
                for k = 1:dof
                    AoplusJ[(i-1)*stageNum+l,(j-1)*stageNum+k] =btab.a[i,j]*J[l,k]
                end
            end
        end
    end

    AoplusI = [btab.a[i,j]*I for i=1:stageNum, j=1:stageNum]
    AoplusI2=zeros(stageNum*dof,stageNum*dof)
    for i=1:stageNum
        for j = 1:stageNum
            for l = 1:dof
                for k = 1:dof
                    AoplusI2[(i-1)*stageNum+l,(j-1)*stageNum+k] =btab.a[i,j]*I[l,k]
                end
            end
        end
    end

    #IoplusM = [I[i,j]*M for i=1:stageNum, j=1:stageNum]
    IoplusM=zeros(stageNum*dof,stageNum*dof)
    for i=1:stageNum
        for j = 1:stageNum
            for l = 1:dof
                for k = 1:dof
                    IoplusM[(i-1)*stageNum+l,(j-1)*stageNum+k] =I[i,j]*M[l,k]
                end
            end
        end
    end

    G =  IoplusM-h*AoplusJ
    Ginv = inv(G)

    # Use Netwon interation
    #
    #   ⃗z^(k+1) = ⃗z ^ (k) - Δ⃗z^(k)
    #TODO: use the transformation T^(-1)A^(-1)T = Λ, W^k = (T^-1 ⊕ I)Z^k
    ## initial variables iteration
    #TODO: use better initial values for zpre
    #w = hnew/hpre
    #zpre = q(w)+ypre-y
    zpre = zeros(dof)
    Δzpre
    κ

    iterate = true
    count = 0
    while iterate
        Δz = reshape(Ginv*[(-zpre + h*AoplusI2*F(f,z,y,t,c,h))...],dof,stageNum)
        z = zpre + Δz

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

    d = b*inv(A)
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
    roots = Array(Float64, stageNum - 1)
    append!(roots, [1 for i= 1:stageNum])
    poly = Polynomials.poly([roots;])
    for i = 1 : stageNum-1
        poly = Polynomials.polyder(poly)
    end
    C = Polynomials.roots(poly)

    ################# Calculate b_i #################

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

    return TableauRKImplicit(order, A, B, C)
end

" Calculates the array of derivative values between t and tnext"
function F(f,z,y,t,c,h)
    return [f(t+c[i]*h, y+z[i]) for i=1:length(z)]
end
