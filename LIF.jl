module LIF

export LIF_spike

using Random

"""
LIF_spike(;V0 =-60e-3, dt =  1e-3, t_max = 400e-3,
              tau= 20e-3, el =-60e-3, i_mean=  25e-3,
              vr =-70e-3, vth=-50e-3, n_amp = 1)

This function implements the Leaky Integrate-and-Fire neuron model.
It requires no input but can optionally receive any of
the model parameters as input to change form the default.
It returns 3 outputs: an array with the time course of the
simulation, an array with the voltage over that time
course, and an array with the spike times. Default parameters
will implement the model with noisy input current whose amplitude
can be modified by modifying the n_amp optional parameter.

# OPTIONAL PARAMETERS
- V0        initial voltage
- dt        time step
- t_max     end time of simulation
- tau       integration constant
- el        leak
- i_mean    mean input current
- vr        resting membrane potential
- vth       threshold voltage for spiking
- n_amp     amplitude of noise

# RETURNS
- time      vector representing the time axis
- v         vector, same size as time, representing voltage
- spike_times   vector of spike times


"""
function LIF_spike(;V0 =-60e-3, dt =  1e-3, t_max = 400e-3,
              tau= 20e-3, el =-60e-3, i_mean=  25e-3,
              vr =-70e-3, vth=-50e-3, vplus= 10e-3, n_amp = 1)

    time = collect(0.0:dt:t_max)
    v = zeros(length(time))
    v[1] = V0
    spike_times = []
    for j in 2:length(time)
        if v[j-1] >= vth
            v[j-1] = vplus
            v[j] = vr
            push!(spike_times,time[j-1])
        else
            i=i_mean*n_amp*rand()
            v[j] = v[j-1] + (el - v[j-1] + i)*dt/tau

        end
    end
    return time, v, spike_times
end #end function

end #end module
