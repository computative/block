import Statistics;
# Returns an estimate for the sample variance of the mean which corrects for correlation in the data
function block(x)
    n = length(x)
    d = log2(n)
    if d%1 != 0
        d = Int(floor(d))
        println("Length of data = $(n) is not a power of 2, truncating to $(Int(round((2^d),digits=0))) elements")
        n = 2^d
        x = x[1:n]
    else
        d = Int(d)
    end
    
    s, gamma = zeros(d), zeros(d)
    mu = Statistics.mean(x)

    # Estimate the auto-covariance and variances for each blocking transformation
    for i in 1:d
        # n changes in length
        n = length(x)

        # Estimate autocovariance of x
        gamma[i] = (1 / n) * sum((x[1:n-1] .- mu) .* (x[2:n] .- mu))

        # Estimate variance of x
        s[i] = Statistics.var(x, corrected=false)

        # Blocking transformation
        x_1 = x[1:2:end] # Extracting all numbers at odd positions
        x_2 = x[2:2:end] # Numbers at even positions
        
        x = 0.5 .* (x_1 .+ x_2)
    end
    
    # Test observator from theorem (chi^2-distributed)
    factor_1 = (gamma ./ s).^2
    factor_2 = 2 .^[i for i in 1:d]
    
    M = (cumsum( ( factor_1 .* factor_2[end:-1:1] )[end:-1:1]))[end:-1:1]

    # Test percentiles
    q = [6.634897,  9.210340, 11.344867, 13.276704, 15.086272,
        16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
        24.724970, 26.216967, 27.688250, 29.141238, 30.577914,
        31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
        38.932173, 40.289360, 41.638398, 42.979820, 44.314105,
        45.641683, 46.962942, 48.278236, 49.587884, 50.892181]

    # The actual Chi squared test - should we have stopped blocking?
    k_end = 1
    for k in 1:d
        k_end = k
        if (M[k] < q[k])
            break
        end
    end
    if (k_end >= d)
        print("More data is needed!")
    end
    result = s[k_end] / 2^(d-(k_end-1))

    return result
end