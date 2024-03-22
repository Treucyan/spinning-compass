# Let's put all of our codes here - Roy 
# Remember to put documentations on all the functions here. - Roy 
# (delete this message before submission)

"""
    EOM_compass(r, t, B, ω)

Equations of motion for a spinning compass in a periodically changing magnetic field

"""
function EOM_compass(r, t, B, ω)
    (x, v) = r
    dxdt = v
    dvdt = -B*cos(ω*t)*sin(x)
    return [dxdt, dvdt]
end    