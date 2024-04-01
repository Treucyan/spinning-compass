# Let's put all of our codes here - Roy 
# Remember to put documentations on all the functions here. - Roy 
# (delete this message before submission)


function EOM_compass(r::Array, t::Float64, B::Float64, ω::Float64)
    """
    EOM_compass(r, t, B, ω)

    #Description
    Equations of motion for a spinning compass in a periodically changing magnetic field

    # Args
        r (Array): Array containing the compass' state at time t. `r` must have the form: `r = [x, v]`
        t (Float64): time t.
        B (Float64): non-dimensional magnetic field amplitude
        ω (Float64): driving frequency

    # Returns
        [dxdt, dvdt] (Array): system's velocity state at t.
    """
    
    (x, v) = r
    dxdt = v
    dvdt = -B*cos(ω*t)*sin(x)
    return [dxdt, dvdt]
end    



# general Runge Kutta algorithm 
# (Sugested Edit: generalize (xpoints, ypoints) to make the function applicable for multidimensional systems)
function RK4(eom_func::Function, time_param::Array, r::Array)
    """
    RK4(f, tpoints, r)

    # Description
    Solve a system of two first-order ODEs using the fourth-order Runge-Kutta method.

    ## Args 
        eom_func (Function): function corresponding to the equations of motion to be integrated.
            The function must follow the format: `eom_func(r, t)`.

        time_param (Array): initial and final time, and number of steps of the integration. `time_param`
            must follow the format: `time_param = [t_initial, t_final, Npoints]`

        r (Array): initial state of the system. r must follow the same ordering of variables used in `eom_func`.

    ## Returns
        (xpoints, ypoints): Arrays containing the IVP solution of the EOM

    """
    #initializing time parameters
    (t_initial, t_final, Npoints) = time_param
    h = abs(t_final - t_initial) / Npoints
    tpoints = range(t_initial, t_final, Npoints)

    #initializing dynamics array 
    xpoints = zeros(Npoints)
    ypoints = zeros(Npoints)


    #runge kutta algorithm
    for (i, t) in enumerate(tpoints)
        xpoints[i] = r[1]
        ypoints[i] = r[2]
        k1 = h*eom_func(r, t)
        k2 = h*eom_func(r .+ 0.5*k1, t + 0.5*h)
        k3 = h*eom_func(r .+ 0.5*k2, t + 0.5*h)
        k4 = h*eom_func(r .+ k3, t + h)
        r = r .+ (k1 + 2*k2 + 2*k3 + k4) / 6
    end
    return xpoints, ypoints
end