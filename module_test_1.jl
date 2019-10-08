# Module Test 1
# Colton Casto

using PyPlot;
using JLD;
using Statistics

# part a: import the function LIF_spike()
push!(LOAD_PATH, pwd()); import LIF: LIF_spike

# part b: run the function with default values, plot voltage, print spikes
times, v, spike_times = LIF_spike();
plot(times, v)
title("Voltage trace as a function of time")
xlabel("Time (seconds)")
ylabel("Voltage (volts)")
println(spike_times)

# part c: overlay different simulation trajectories with various mean current
times, v_1, spike_times_1 = LIF_spike(i_mean = 0.01);
times, v_2, spike_times_2 = LIF_spike(i_mean = 0.025);
times, v_3, spike_times_3 = LIF_spike(i_mean = 0.05);

clf()
plot(times, v_1)
plot(times, v_2)
plot(times, v_3)
title("Voltage trace with different mean input current values")
xlabel("Time (seconds)")
ylabel("Voltage (volts)")
legend(["mean current = 0.01", "mean current = 0.025", "mean current = 0.05"])


# part d: write a function that runs the simulation N times and calcs mean ISI
"""
avg_ISI(n; i_mean) - This function runs lIF_spike() n times and calculates the
average interspike interval (ISI) for each of the simulations. This function
assumes that if there is only one spike (or fewer) in the time interval that
there is not enough information to calculate an average duration between two
action and will place an arbitrarily large duration for those instances.

Args:
- n (int): number of simulations to run
- i_mean (float): OPTIONAL average input current. Default 0.025.

Returns:
- n_averages (array): array with average ISI for each of the n simulations
"""
function avg_ISI(n; i_mean = 0.025)
    n_averages = zeros(n)
    for i = 1:n
        times, v, spike_times = LIF_spike(i_mean = i_mean)

        # calculate the durations between every spike
        if length(spike_times) > 1 # case with more than one spike
            spike_duration = spike_times[2:length(spike_times)] -
            spike_times[1:length(spike_times)-1]
            n_averages[i] = mean(spike_duration)
        else # case with one spike or no spikes
            n_averages[i] = times[end] + 1 # arbitrarily large number
        end
    end
    return n_averages
end

n = 100
n_averages = avg_ISI(n)
clf()
hist(n_averages, bins = 20)
title("Histogram of average interspike intervals (ISI) for 100 simulations and
an average current input of 0.025")
xlabel("duration")
ylabel("frequency")

# part e: same as (d) but i_mean = 0.045 and 0.05, plot histogram of results
n_averages_45 = avg_ISI(n, i_mean = 0.045)
n_averages_5 = avg_ISI(n, i_mean = 0.05)
clf()
hist(n_averages_45, bins = 20)
hist(n_averages_5, bins = 20)
title("Histogram of average interspike intervals (ISI) for 100 simulations with
an average current input of 0.045 and 0.05")
xlabel("duration")
ylabel("frequency")
legend(["mean current = 0.045", "mean current = 0.05"])

# part f: firing rate vs. input current
"""
f_i_curve(v) - This function calculates the firing rate given different average current input values. It also uses the avg_ISI function to help smooth the curve by running multiple simulations for each mean current.

Args:
- v (array): mean input current values to calculate the firing rate for

Returns:
- firing_rate (array): firing rate for each of the inpu current values given as input
"""
function f_i_curve_new(v)
    firing_rate = zeros(length(v))
    for i = 1:length(v)
        times, v, spike_times = LIF_spike(i_mean = v[i])
        firing_rate[i] = (length(spike_times)/times[end])
    end
    return firing_rate
end

"""
f_i_curve(v) - This function calculates the firing rate given different average current input values. It also uses the avg_ISI function to help smooth the curve by running multiple simulations for each mean current.

Args:
- v (array): mean input current values to calculate the firing rate for

Returns:
- firing_rate (array): firing rate for each of the inpu current values given as input
"""
function f_i_curve_original(v)
    n = 100
    length_t = 0.4
    n_averages = zeros(n)
    firing_rate = zeros(length(v))
    for i = 1:length(v)
        n_averages = avg_ISI(n, i_mean = v[i])
        avg_duration = mean(n_averages)
        firing_rate[i] = length_t/avg_duration
    end
    return firing_rate
end


v = 0:0.001:0.12
firing_rate = f_i_curve_original(v)
clf()
plot(v, firing_rate)
title("Firing Rate vs. Input Current")
xlabel("Input Current")
ylabel("Firing Rate")

# part g: traces overlayed 15ms before spike with optional inputs to LIF_spike
clf()
length_t = 15
times, v, spike_times = LIF_spike(t_max = 0.9, dt = 0.001) # specified in q
spike_instances = zeros(length_t, length(spike_times)) # eg 15x5 matrix
indices = zeros(Int64, length(spike_times))
for i = 1:length(spike_times)
    indices[i] = findall(times .== spike_times[i])[1]
    index = indices[i]
    index_start = index-length_t
    spike_instances[:,i] = v[index_start:index-1]
    plot(-0.015:0.001:-0.001, spike_instances[:,i], alpha = 0.25, color = "k")
end

# find the average trace
mean_trace = zeros(length_t)
for i = 1:length_t
    mean_trace[i] = mean(spike_instances[i,:])
end

plot(-0.015:0.001:-0.001, mean_trace, linewidth = 2, color = "r", alpha = 1)
title("Plot of voltage traces 15ms before spike
(red line = average across spikes)")
xlabel("time in seconds (0 = spike)")
ylabel("voltage (volts)")

# part h: standard error
"""
serror(n; i_mean) - Calculates the mean, variance, and standard error of the
output of avg_ISI().

Args:
- n (int): number of simulations
- i_mean (float): OPTIONAL average input current. Default 0.1.

Returns:
- sem (float): standard error of the mean
"""
function serror(n; i_mean = 0.1)
    n_averages = avg_ISI(n, i_mean = i_mean)
    mu = mean(n_averages) # mean

    # variance calculation
    var = 0
    for i = 1:n
        var = var + ((n_averages[i] - mu) * (n_averages[i] - mu))
        var = var / (n - 1)
    end

    # sem calculation
    return sqrt(var / n)
end

serror(100, i_mean = 0.3)

# part i: standard error plotter
"""
serror_plotter(n; i_mean) - plots standard error of the mean as a function of the number of simulations

Args:
- n (int): number of simulations to run when serror() is called
- i_mean (float): OPTIONAL average input current

Returns:
- None: produces the N vs. SEM graph
"""
function serror_plotter(n; i_mean = 0.1)
    clf()
    sem = zeros(length(n))
    for i = 1:length(n)
        sem[i] = serror(n[i], i_mean = i_mean)
    end

    plot(n, sem)
    title("Standard error as a function of the number of simulations")
    xlabel("N")
    ylabel("Standard error")
end

N = [20, 40, 80, 160, 320]
serror_plotter(N)
