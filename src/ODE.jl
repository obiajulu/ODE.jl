isdefined(Base, :__precompile__) && __precompile__()
# Ordinary Differential Equation Solvers

module ODE

using Polynomials
using Compat
using Iterators

## minimal function export list
# adaptive non-stiff:
export ode23, ode45, ode78
# non-adaptive non-stiff:
export ode4, ode4ms
# adaptive stiff:
export ode23s
# non-adaptive stiff:
export ode4s

import Base.convert, Base.show
import Base: start, next, done, call, collect

## complete function export list: see runtests.jl

<<<<<<< HEAD
    # initialization
    t = tspan[1]

    tfinal = tspan[end]

    h = initstep
    if h == 0.
        # initial guess at a step size
        h, tdir, F0 = hinit(F, y0, t, tfinal, 3, reltol, abstol)
    else
        tdir = sign(tfinal - t)
        F0 = F(t,y0)
    end
    h = tdir * min(abs(h), maxstep)

    y = y0
    tout = Array(typeof(t), 1)
    tout[1] = t         # first output time
    yout = Array(typeof(y0), 1)
    yout[1] = deepcopy(y)         # first output solution


    J = jac(t,y)    # get Jacobian of F wrt y

    while abs(t - tfinal) > 0 && minstep < abs(h)
        if abs(t-tfinal) < abs(h)
            h = tfinal - t
        end

        if size(J,1) == 1
            W = one(J) - h*d*J
        else
            # note: if there is a mass matrix M on the lhs of the ODE, i.e.,
            #   M * dy/dt = F(t,y)
            # we can simply replace eye(J) by M in the following expression
            # (see Sec. 5 in [SR97])

            W = lufact( eye(J) - h*d*J )
        end

        # approximate time-derivative of F
        T = h*d*(F(t + h/100, y) - F0)/(h/100)

        # modified Rosenbrock formula
        k1 = W\(F0 + T)
        F1 = F(t + 0.5*h, y + 0.5*h*k1)
        k2 = W\(F1 - k1) + k1
        ynew = y + h*k2
        F2 = F(t + h, ynew)
        k3 = W\(F2 - e32*(k2 - F1) - 2*(k1 - F0) + T )

        err = (abs(h)/6)*norm(k1 - 2*k2 + k3) # error estimate
        delta = max(reltol*max(norm(y),norm(ynew)), abstol) # allowable error

        # check if new solution is acceptable
        if  err <= delta

            if points==:specified || points==:all
                # only points in tspan are requested
                # -> find relevant points in (t,t+h]
                for toi in tspan[(tspan.>t) & (tspan.<=t+h)]
                    # rescale to (0,1]
                    s = (toi-t)/h

                    # use interpolation formula to get solutions at t=toi
                    push!(tout, toi)
                    push!(yout, y + h*( k1*s*(1-s)/(1-2*d) + k2*s*(s-2*d)/(1-2*d)))
                end
            end
            if (points==:all) && (tout[end]!=t+h)
                # add the intermediate points
                push!(tout, t + h)
                push!(yout, ynew)
            end

            # update solution
            t = t + h
            y = ynew

            F0 = F2         # use FSAL property
            J = jac(t,y)    # get Jacobian of F wrt y
                            # for new solution
        end

        # update of the step size
        h = tdir*min( maxstep, abs(h)*0.8*(delta/err)^(1/3) )
    end

    return tout, yout
end


#ODEROSENBROCK Solve stiff differential equations, Rosenbrock method
#   with provided coefficients.
function oderosenbrock(F, x0, tspan, gamma, a, b, c; jacobian=nothing)

    if typeof(jacobian) == Function
        G = jacobian
    else
        G = (t, x)->fdjacobian(F, x, t)
    end

    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    solstep = 1
    while solstep < length(tspan)
        ts = tspan[solstep]
        hs = h[solstep]
        xs = x[solstep]
        dFdx = G(ts, xs)
        # FIXME
        if size(dFdx,1) == 1
            jac = 1/gamma/hs - dFdx[1]
        else
            jac = eye(dFdx)/gamma/hs - dFdx
        end

        g = Array(typeof(x0), size(a,1))
        g[1] = (jac \ F(ts + b[1]*hs, xs))
        x[solstep+1] = x[solstep] + b[1]*g[1]

        for i = 2:size(a,1)
            dx = zero(x0)
            dF = zero(x0/hs)
            for j = 1:i-1
                dx += a[i,j]*g[j]
                dF += c[i,j]*g[j]
            end
            g[i] = (jac \ (F(ts + b[i]*hs, xs + dx) + dF/hs))
            x[solstep+1] += b[i]*g[i]
        end
        solstep += 1
    end
    return vcat(tspan), x
end


# Kaps-Rentrop coefficients
const kr4_coefficients = (0.231,
                          [0              0             0 0
                           2              0             0 0
                           4.452470820736 4.16352878860 0 0
                           4.452470820736 4.16352878860 0 0],
                          [3.95750374663  4.62489238836 0.617477263873 1.28261294568],
                          [ 0               0                0        0
                           -5.07167533877   0                0        0
                            6.02015272865   0.1597500684673  0        0
                           -1.856343618677 -8.50538085819   -2.08407513602 0],)

ode4s_kr(F, x0, tspan; jacobian=nothing) = oderosenbrock(F, x0, tspan, kr4_coefficients...; jacobian=jacobian)

# Shampine coefficients
const s4_coefficients = (0.5,
                         [ 0    0    0 0
                           2    0    0 0
                          48/25 6/25 0 0
                          48/25 6/25 0 0],
                         [19/9 1/2 25/108 125/108],
                         [   0       0      0   0
                            -8       0      0   0
                           372/25   12/5    0   0
                          -112/125 -54/125 -2/5 0],)

ode4s_s(F, x0, tspan; jacobian=nothing) = oderosenbrock(F, x0, tspan, s4_coefficients...; jacobian=jacobian)
=======
# basic type definitions
include("types.jl")
include("helpers.jl")
>>>>>>> pwl/master

# dense output wrapper
include("dense.jl")

# particular solvers
include("ode23s.jl")
include("rk.jl")
# include("multistep.jl")

include("iterators.jl")
include("interfaces.jl")

end # module ODE
