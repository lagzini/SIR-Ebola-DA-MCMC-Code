function data_set = Gillespie_SIR2(N, beta, gamma, i0, tmax)
    n = N;
    s0 = n - i0;
    time = 0;
    S = s0;
    I = i0;
    R = n - s0 - i0;

    data_set = [time, I, S, R];

    while time(end) < tmax
        rate1 = beta * S(end) * I(end) / n;
        rate2 = gamma * I(end);
        total_rate = rate1 + rate2;

        tau = exprnd(1 / total_rate);

        if rand < rate1 / total_rate
            S = [S, S(end) - 1];
            I = [I, I(end) + 1];
            R = [R, R(end)];
        else
            I = [I, I(end) - 1];
            R = [R, R(end) + 1];
            S = [S, S(end)];
        end

        time = [time, time(end) + tau];
        data_set = [data_set; time(end), I(end), S(end), R(end)];
    end
end
