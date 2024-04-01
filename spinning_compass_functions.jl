# Let's put all of our codes here - Roy 
# Remember to put documentations on all the functions here. - Roy 
# (delete this message before submission)




#module for the equations of motions and rk4 algorithm
module Spin_compass


#equations of motions of the spinning compass
"
EOM_compass(r, t, B, ω)

#Description
Equations of motion for a spinning compass in a periodically changing magnetic field

# Args
    r (Array): Array containing the compass' state at time t. `r` must have the form: `r = [x, v]`
    t (Float64): time t.
    B (Float64): non-dimensional magnetic field amplitude
    ω (Float64): non-dimensional driving frequency

# Returns
    [dxdt, dvdt] (Array): system's velocity state at t.
"
function EOM_compass(r::Array, t::Float64, B::Float64, ω::Float64)
    
    
    (x, v) = r
    dxdt = v
    dvdt = -B*cos(ω*t)*sin(x)
    return [dxdt, dvdt]
end    



# general Runge Kutta algorithm 
# (Sugested Edit: generalize (xpoints, ypoints) to make the function applicable for multidimensional systems)


"
RK4(f, time_param, r)

# Description
Solve a system of two first-order ODEs using the fourth-order Runge-Kutta method.

## Args 
    eom_func (Function): function corresponding to the equations of motion to be integrated.
        The function must follow the format: `eom_func(r, t)`.

    time_param (Tuple): initial and final time, and number of steps of the integration. `time_param`
        must follow the format: `time_param = [t_initial, t_final, Npoints]`

    r (Array): initial state of the system. r must follow the same ordering of variables used in `eom_func`.

## Returns
    (xpoints, ypoints): Arrays containing the IVP solution of the EOM

"
function RK4(eom_func::Function, time_param::Tuple{Float64, Float64, Int64}, r::Array)
    
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
    return tpoints, xpoints, ypoints
end


"
module_initializer()

# Descrition
    Prints Hello World!
"
function module_initializer()
    println("Hello World")
end

end


module Chaos_checking

using FFTW

#spectral entropy
"spectral_entropy(observable)

# Description
Calculates the spectral entropy of a single observable. Note that the function has a dependency on the FFTW package
as it uses the fft function to do a fast fourier transform.
"
function spectral_entropy(observable::Array)
    fourier_observable = fftshift(fft(observable))
    power_spectrum = abs.(fourier_observable).^2
    normalized_power = (power_spectrum./ sum(power_spectrum)).*2

    #calculating spectral entropy
    H_spectral = 0
    for probability in normalized_power
        if probability > 1e-5
            H_spectral += probability * log2(probability)
        else
            H_spectral += 0
        end
    end
    H_spectral /= -2* log2(length(normalized_power))
    return H_spectral
end





end