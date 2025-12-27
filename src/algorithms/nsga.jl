function obj_NSGA(PROBLEM::String)

    if PROBLEM == "F1" || PROBLEM == "F7" || PROBLEM == "F8"
        n = 10
        m = 2
    elseif PROBLEM == "F6"
        n = 10
        m = 3
    elseif PROBLEM == "F2" || PROBLEM == "F3" || PROBLEM == "F4" || 
           PROBLEM == "F5" || PROBLEM == "F9"
        n = 30
        m = 2
    elseif PROBLEM == "F6"
        n = 10
        m = 3
    elseif PROBLEM == "FDS"
        n = 200
        m = 3
    elseif PROBLEM == "JOS1"
        n = 100
        m = 2
    elseif PROBLEM == "QV1"
        n = 100
        m = 2
    elseif PROBLEM == "SLCDT2"
        n = 100
        m = 3
    elseif PROBLEM == "ZDT1" || PROBLEM == "ZDT2" || PROBLEM == "ZDT3"
        n = 30
        m = 2
    elseif PROBLEM == "ZDT4" || PROBLEM == "ZDT6"
        n = 10
        m = 2
    elseif PROBLEM == "ZLT1"
        n = 100
        m = 5
    else
        error("Problem not yet implemented for NSGAII: $PROBLEM")
    end

    return function (x)
        F = zeros(m)
        for i in 1:m
            F[i] = evalf(i, x, n)
        end
        return F, [0.0], [0.0]
    end

end
