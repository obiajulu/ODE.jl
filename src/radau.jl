#=
    (1) done
    (2) trialstep!()
    (3) errorcontol!()
    (4) odercontrol!()
    (5) rollback!()
    (6) status()
=#

function radau(f, y0, tspan, order ::Integer = 5)
    #= setup
        state.t, state.f(t)
    =#

    while !done()
        stats = trialstep!(odep, st)
        err, stats, st = errorcontrol!(odep, st) # (1) adaptive stepsize (2) error
        if err < 1
            stats, st = ordercontrol!()
            accept = true
        else
            rollback!()
        end

        return = status()
    end
end

function done(arg..)
