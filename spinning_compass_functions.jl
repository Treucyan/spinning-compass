# Let's put all of our codes here - Roy 
# Remember to put documentations on all the functions here. - Roy 
# (delete this message before submission)




# --------------------------------------------------
#  Spin_compass
# --------------------------------------------------
module Spin_compass

"""
EOM_compass(r, t, B, ω)

# Description
Equations of motion for a spinning compass in a periodically changing magnetic field

# Args
    r (Array): Array containing the compass' state at time t. `r` must have the form: `r = [x, v]`
    t (Float64): time t.
    B (Float64): non-dimensional magnetic field amplitude
    ω (Float64): non-dimensional driving frequency

# Returns
    [dxdt, dvdt] (Array): system's velocity state at t.
"""
function EOM_compass(r::Array, t::Float64, B::Float64, ω::Float64)
    (x, v) = r
    dxdt = v
    dvdt = -B*cos(ω*t)*sin(x)
    return [dxdt, dvdt]
end    


"""
EOM_compass_unitless(r, t, λ)

# Description 
Unitless equations of motions of the periodically-driven spinning compass. 

The time variable is expressed in terms of the driving period `T = ωt``. For the simulation time to be commensurate with the driving period, set T_f = `2πn`

# Args
    r (Array): Array containing the compass' state at time t. `r` must have the form: `r = [x, v]`
    t (Float64): time t
    λ (Float64): unitless amplitude

# Returns
    [dxdt, dvdt] (Array): system's velocity state at t.
"""
function EOM_compass_unitless(r::Array, t::Float64, λ::Float64)
    (x, v) = r
    dxdT = v
    dvdT = -0.5 * (λ^2) * cos(t) * sin(x)
    return [dxdT, dvdT]
end


"""
    RK4(eom_func, t_param, r)

Solve systems of first-order ordinary differential equations using the fourth-order Runge-Kutta method.

# Args
- eom_func (Function): function corresponding to the equations of motion to be integrated. The function must follow the format: `eom_func(r, t)`.

- time_param (Tuple): initial and final time, and number of steps of the integration. `time_param` must follow the format: `time_param = [t_initial, t_final, Npoints]`

- r (Array): initial state of the system. r must follow the same ordering of variables used in `eom_func`.

# Returns
- (tpoints, rpoints): Arrays containing the IVP solution of the EOM

# Notes
- Each row in rpoints contains the values of a particular dependent variable. To extract these values, one can either do
`xpoints = rpoints[1, :], ypoints = rpoints[2, :], ...` or `xpoints, ypoints, ... = eachrow(rpoints).`
"""
function RK4(eom_func::Function, t_param::Tuple{Float64, Float64, Int64}, r::Vector{Float64})

    # Initialize time array
    (t_initial, t_final, Npoints) = t_param
    h = abs(t_final - t_initial) / Npoints
    tpoints = range(t_initial, t_final, Npoints)

    # Initialize trajectory array
    rpoints = zeros(length(r), Npoints)

    # Runge-Kutta algorithm
    for (i, t) in enumerate(tpoints)
        rpoints[:, i] = r
        k1 = h*eom_func(r, t)
        k2 = h*eom_func(r .+ (0.5 .* k1), t + 0.5*h)
        k3 = h*eom_func(r .+ (0.5 .* k2), t + 0.5*h)
        k4 = h*eom_func(r .+ k3, t + h)
        r = r .+ (k1 + 2 .* k2 + 2 .* k3 + k4) / 6
    end
    return tpoints, rpoints
end


"""
module_initializer()

# Descrition
    Prints Hello World!
"""
function module_initializer()
    println("Hello World")
end

end





# --------------------------------------------------
#  Chaos_checking
# --------------------------------------------------
module Chaos_checking
using FFTW


"""
    spectral_entropy(observable)

# Description
Calculates the spectral entropy of a single observable. Note that the function has a dependency on the FFTW package
as it uses the fft function to do a fast fourier transform.
"""
function spectral_entropy(x, cutoff::Float64 = 1e-10)
    fourier_x = fftshift(fft(x))
    power_spectrum = abs2.(fourier_x)  # abs2 slightly faster than abs()^2
    power_spectrum ./= sum(power_spectrum)  # normalization

    spectral = filter(x -> x > cutoff, power_spectrum)  # remove items below cutoff
    spectral .*= log2.(spectral)
    spectral /= -1 * log2(length(power_spectrum))

    return sum(spectral)
end


"""
spectral_entropy(observable)

# Description
Calculates the spectral entropy of a single observable. Note that the function has a dependency on the FFTW package
as it uses the fft function to do a fast fourier transform.
"""
function spectral_entropy_old(observable::Array)
    fourier_observable = fftshift(fft(observable))
    power_spectrum = abs.(fourier_observable).^2 
    normalized_power = power_spectrum./ sum(power_spectrum)

    #calculating spectral entropy
    H_spectral = 0
    @simd for probability in normalized_power
        if probability >= 1e-10
            H_spectral += probability * log2(probability)
        else
            H_spectral += 0
        end
    end
    H_spectral /= -1 * log2( length(normalized_power))
    return H_spectral
end


"""
stroboscope_dynamics(observables, T_array, timestep)

# Description
Constructs the stroboscopic dynamics of a given set of observables, `observables = (xpoints, vpoints)`.

## Args
    xpoints (Array): full dynamics of the position variable
    
    vpoints (Array): full dynamics of the velocity variable

    time_param (Tuple): time parameters for constructing the time array. Use the following
        format: time_param = [t_initial, t_final, Nsteps]

## Returns 
    x_strobe, v_strobe (Vector): stroboscopic dynamics of xpoints and vpoints
"""
function stroboscope_dynamics(xpoints, vpoints, time_param::Tuple{Float64, Float64, Int64})
    (t_initial, t_final, Nsteps) = time_param
    timestep = abs(t_final - t_initial) / Nsteps
    T_array = range(t_initial, t_final, Nsteps)

    x_strobe, v_strobe = [], []
    @simd for i in 1:Nsteps
        time_error = T_array[i] % 2π
        if time_error <= timestep
            append!(x_strobe, xpoints[i])
            append!(v_strobe, vpoints[i])
        end
    end
    return x_strobe, v_strobe
end




end



# --------------------------------------------------
#  Phase_diagram
# --------------------------------------------------
module Phase_diagram
using DelimitedFiles
using FFTW
using ..Spin_compass
using ..Chaos_checking


"""
lambda_entropy_linear_scan(time_param, scan_param)

# Description
Generates a linear scan of lambda for the unitless spin compass, calculating the 
spectral entropy at each value of lambda.

## Args
    time_param (Tuple{Float64, Float64, Int64}): initial and final time, and number of steps of the integration. 
        `time_param` must follow the format: `time_param = [t_initial, t_final, Npoints]`

    scan_param (Tuple{Float64, Float64, Int64}): range of lambda to be considered and the resolution of the scan.
        `scan_param` must follow the format: 'scan_param = [lambda_initial, lambda_final, resolution]'

## Returns
    spec_entropy_array (Array): spectral entropy as a function of lambda for the given range set by `scan_param`
"""
function lambda_entropy_linear_scan(time_param::Tuple{Float64, Float64, Int64}, scan_param::Tuple{Float64, Float64, Int64})
    #initializing constants
    (lambda_initial, lambda_final, resolution) = scan_param

    #scan Arrays
    lambda_sample_array = range(lambda_initial, lambda_final, resolution)
    spec_entropy_array = zeros(resolution)

    println("Scan starting...")
    for i in 1:resolution
        #initializing dynamics inputs
        λ = lambda_sample_array[i]
        x0, v0 = 1.0, 0.0
        r = [x0, v0]
        f(r, t) = Spin_compass.EOM_compass_unitless(r, t, λ)

        (tpoints, rpoints) = Spin_compass.RK4(f, time_param, r)
        xpoints = rpoints[1, :]

        #bounding phi
        cartesian_proj_x = cos.(xpoints)
        cartesian_proj_y = sin.(xpoints)
        xpoints = atan.(cartesian_proj_y, cartesian_proj_x) # we do this to force phi to be periodic around -pi and pi

        entropy = Chaos_checking.spectral_entropy(xpoints)

        spec_entropy_array[i] = entropy

        if i % 10 == 0
            println("Number of Rows Done: ", i)
            println("Number of Rows Remaining: ", resolution - i)
        end
    end

    println("Scan completed!")

    return spec_entropy_array
end


"""
lambda_linear_scan_saver(time_param, scan_param, save_filename)

# Description
Generates and saves the spectral entropy array generated from the linear scan of lambda into a txt file. Note that
the function uses the function `writedlm()` under the package `DelimitedFiles`.

## Args
    time_param (Tuple{Float64, Float64, Int64}): initial and final time, and number of steps of the integration. 
        `time_param` must follow the format: `time_param = [t_initial, t_final, Npoints]`

    scan_param (Tuple{Float64, Float64, Int64}): range of lambda to be considered and the resolution of the scan.
        `scan_param` must follow the format: 'scan_param = [lambda_initial, lambda_final, resolution]'

    save_filename (String): filename of the txt where the spectral entroy array will be saved

## Returns
    A txt file containing the spectral entropy array for the given range of lambda set by `scan_param`

"""
function lambda_linear_scan_saver(time_param::Tuple{Float64, Float64, Int64}, scan_param::Tuple{Float64, Float64, Int64}, save_filename::String)
    spec_entropy_array = Phase_diagram.lambda_entropy_linear_scan(time_param, scan_param)
    writedlm(save_filename, spec_entropy_array)
end


"""
normalized_power_linear_scan(time_param, scan_param)

# Description
Generates a linear scan of lambda to obtain the behavior of the normalized power spectrum
of the spin compass as a function of lambda.

## Args
    time_param (Tuple{Float64, Float64, Int64}): initial and final time, and number of steps of the integration. 
        `time_param` must follow the format: `time_param = [t_initial, t_final, Npoints]`

    scan_param (Tuple{Float64, Float64, Int64}): range of lambda to be considered and the resolution of the scan.
        `scan_param` must follow the format: 'scan_param = [lambda_initial, lambda_final, resolution]'

## Returns
    freq_spectrum (Array): angular frequency axis of the power spectrum
    normalized_power_fullarray (Array): 2x2 matrix of the normalized power spectrum as a function of lambda

"""
function normalized_power_linear_scan(time_param::Tuple{Float64, Float64, Int64}, scan_param::Tuple{Float64, Float64, Int64})
    x0, v0 = 1.0, 0.0
    r = [x0, v0]

    #initializing constants
    (t_initial, t_final, Nsteps) = time_param
    (lambda_initial, lambda_final, resolution) = scan_param
    sampling_rate = Nsteps / abs(t_final - t_initial)

    #scan Arrays
    lambda_sample_array = range(lambda_initial, lambda_final, resolution)
    normalized_power_fullarray = zeros((Nsteps, resolution))
    freq_spectrum = fftshift(fftfreq(Nsteps, sampling_rate)) * 2π

    println("Scan starting...")
    for i in 1:resolution
        #initializing dynamics inputs
        λ = lambda_sample_array[i]
        f(r, t) = Spin_compass.EOM_compass_unitless(r, t, λ)

        (tpoints, rpoints) = Spin_compass.RK4(f, time_param, r)
        xpoints = rpoints[1, :]

        #calculating the frequency spectrum
        #bounding phi
        cartesian_proj_x = cos.(xpoints)
        cartesian_proj_y = sin.(xpoints)
        xpoints = atan.(cartesian_proj_y, cartesian_proj_x) # we do this to force phi to be periodic around -pi and pi

        fourier_xpoints = fftshift(fft(xpoints))
        power_spectrum = abs2.(fourier_xpoints)
        normalized_power_fullarray[:, i] = power_spectrum ./ sum(power_spectrum)

        #for code progress tracking
        if i % 10 == 0
            println("Number of Rows Done: ", i)
            println("Number of Rows Remaining: ", resolution - i)
        end
    end


    println("Scan completed!")

    return freq_spectrum, normalized_power_fullarray
end


"""
normalized_power_scan_saver(time_param, scan_param, freq_filename, normalized_power_filename)

# Description
Generates and saves the normalized power spectrum array generated from the linear scan of lambda into a txt file. Note that
the function uses the function `writedlm()` under the package `DelimitedFiles`.

## Args
    time_param (Tuple{Float64, Float64, Int64}): initial and final time, and number of steps of the integration. 
        `time_param` must follow the format: `time_param = [t_initial, t_final, Npoints]`   

    scan_param (Tuple{Float64, Float64, Int64}): range of lambda to be considered and the resolution of the scan.
        `scan_param` must follow the format: 'scan_param = [lambda_initial, lambda_final, resolution]'

    freq_filename (String): filename of the txt file where the frequency spectrum array will be saved

    normalized_power_filename (String): filename of the txt file where the power spectrum array will be saved

## Returns
    A txt file containing the power spectrum array for the given range of lambda set by `scan_param`

"""
function normalized_power_scan_saver(time_param::Tuple{Float64, Float64, Int64}, scan_param::Tuple{Float64, Float64, Int64}, freq_spectrum_filename::String, normalized_power_filename::String)
    freq_spectrum, normalized_power = normalized_power_linear_scan(time_param, scan_param)
    writedlm(freq_spectrum_filename, freq_spectrum)
    writedlm(normalized_power_filename, normalized_power)
end


"""
b_omega_spectral_diagram_scanner(time_param, scan_param)

# Description
Generates a 2D array with size (resolution , resolution) containing the spectral entropy of the spin compass' 
dynamics as a function of both B and ω for a given range set by the `scan_param` input.

## Args
    scan_param (Tuple{Float64, Float64, Float64, Float64, Int64}): range of lambda to be considered and the resolution of 
        the scan. `scan_param` must follow the format: 
        'scan_param = [B_initial, B_final, omega_initial, omega_final, resolution]'

## Returns
    spectral_entropy_diagram (Array): 2D array containing the spectral entropy as a function of B and ω
"""
function b_omega_spectral_diagram_scanner(scan_param::Tuple{Float64, Float64, Float64, Float64, Int64})
    (B_initial, B_final, omega_initial, omega_final, resolution) = scan_param

    #initializing Arrays
    magnetic_field_amp_array = range(B_initial, B_final, resolution)
    driving_freq_array = range(omega_initial, omega_final, resolution)
    spectral_entropy_diagram = zeros((resolution, resolution))

    #initial states
    x0, v0 = 1.0, 0.0
    r = [x0, v0]

    timestep = 0.01

    println("Scan starting...")
    for j in 1:resolution

        ω = driving_freq_array[j]
        t_initial, t_final = 0.0, 100 * 2π / ω
        Nsteps = Int64(floor(abs(t_final - t_initial) / timestep))
        time_param = (t_initial, t_final, Nsteps)

        for i in 1:resolution
            B = magnetic_field_amp_array[i]
            f(r, t) = Spin_compass.EOM_compass(r, t, B, ω)

            (tpoints, rpoints) = Spin_compass.RK4(f, time_param, r)
            xpoints = rpoints[1, :]

            #bounding phi
            cartesian_proj_x = cos.(xpoints)
            cartesian_proj_y = sin.(xpoints)
            xpoints = atan.(cartesian_proj_y, cartesian_proj_x) # we do this to force phi to be periodic around -pi and pi

            spectral_entropy_diagram[i, j] = Chaos_checking.spectral_entropy(xpoints)
        end
        if j % 10 == 0
            println("Number of Columns Done: ", j)
            println("Number of Columns Remaining: ", resolution - j)
        end
    end

    println("Scan completed!")

    return spectral_entropy_diagram
end


"""
b_omega_spectral_diagram_scan_saver()
# Description
Generates and saves the matrix containing the spectral entropy information as a function of B and ω information
a txt file.

## Args
    scan_param (Tuple{Float64, Float64, Float64, Float64, Int64}): range of lambda to be considered and the resolution of 
        the scan. `scan_param` must follow the format: 
        'scan_param = [B_initial, B_final, omega_initial, omega_final, resolution]'

    save_filename (String): filename of the txt file where the spectral entropy diagram will be saved


## Returns
    A txt file containing the spectral entropy diagram as a function of B and ω for the given range set by scan_param.

"""
function b_omega_spectral_diagram_scan_saver(scan_param::Tuple{Float64, Float64, Float64, Float64, Int64}, save_filename::String)
    spectral_entropy_diagram = Phase_diagram.b_omega_spectral_diagram_scanner(scan_param)
    writedlm(save_filename, spectral_entropy_diagram)
end




end