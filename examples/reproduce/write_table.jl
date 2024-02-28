using Query

function save_binned_data(fn, df; bin_pos_emt = false, Δt = 6.0, Δy = 0.333)

    # cleanup
    map!( s -> replace(s, "R" => ""), df.scenario, df.scenario)
    map!( s -> replace(s, "P" => "S"), df.scenario, df.scenario)
    map!( s -> replace(s, "∅" => "0"), df.scenario, df.scenario)

    scenarios = unique(df.scenario)
    sort_scenarios!(scenarios)

    dfd = deepcopy(df)
    function round_dt(x, Δt)
        if isinf(Δt)
            return if isinf(x); Inf64 else 0.0 end
        else
            return Δt * (round( x / Δt + 0.5 ) - 0.5)
        end
    end
    function first_time(row)
        i = argmin([row.A, row.B, row.S])
        return [row.y_A, row.y_B, row.y_S][i]
    end

    dfd[!,:pos_emt] = first_time.(eachrow(dfd))
    dfd[!,:pos_emt] .= round_dt.(dfd.pos_emt, Δy)

    dfd.A .= round_dt.(dfd.A, Δt)
    dfd.B .= round_dt.(dfd.B, Δt)
    dfd.S .= round_dt.(dfd.S, Δt)

    groups = groupby(dfd, [:INM,:A,:B,:S,:pos_emt,:running_mode])
    sat = [nrow(g) for g in groups]
    median(sat)
    mean(sat)
    minimum(sat)
    printstyled("""Number of cases in bins (Δt = $Δt): mean = $(mean(sat)), median = $(median(sat)), minimum = $(minimum(sat))\n""",
                    color = :green)

    mean_or_nan(x) = isempty(x) ? NaN64 : mean(x)
    std_or_nan(x) = isempty(x) ? NaN64 : std(x)

    if bin_pos_emt
        df2 = dfd |>
            @groupby((_.scenario,_.INM,_.running_mode,_.pos_emt,_.A,_.B,_.S)) |>
            @map(  {
                    scenario=key(_)[1],
                    INM=Int64(key(_)[2]),
                    running_mode=key(_)[3],
                    pos_emt=key(_)[4],
                    A=key(_)[5],
                    B=key(_)[6],
                    S=key(_)[7],
                    t_basal_extr_mean = mean_or_nan(t for t in _.t_basal_extr if t >= 0.0),
                    t_basal_extr_std = std_or_nan(t for t in _.t_basal_extr if t >= 0.0),
                    n = length(_),
                    basal_extr_count = count(_.basal_extr),
                    apical_extr_count = count(_.apical_extr)}) |>
            @mutate( basal_extr_ratio = _.basal_extr_count / _.n,
                     apical_extr_ratio = _.apical_extr_count / _.n) |>
            @orderby(_.apical_extr_ratio) |>
            @orderby_descending(_.basal_extr_ratio) |>
            DataFrame


        CSV.write(fn * ".csv", df2; decimal = ',', delim='\t')
        XLSX.writetable(fn * ".xlsx", collect(DataFrames.eachcol(df2)), DataFrames.names(df2), overwrite = true)
    else
        df2 = dfd |>
            @groupby((_.scenario,_.INM,_.running_mode,_.A,_.B,_.S)) |>
            @map(  {
                    scenario=key(_)[1],
                    INM=Int64(key(_)[2]),
                    running_mode=key(_)[3],
                    A=key(_)[4],
                    B=key(_)[5],
                    S=key(_)[6],
                    t_basal_extr_mean = mean_or_nan(t for t in _.t_basal_extr if t >= 0.0),
                    t_basal_extr_std = std_or_nan(t for t in _.t_basal_extr if t >= 0.0),
                    n = length(_),
                    basal_extr_count = count(_.basal_extr),
                    apical_extr_count = count(_.apical_extr)}) |>
            @mutate( basal_extr_ratio = _.basal_extr_count / _.n,
                     apical_extr_ratio = _.apical_extr_count / _.n) |>
            @orderby(_.apical_extr_ratio) |>
            @orderby_descending(_.basal_extr_ratio) |>
            DataFrame

        CSV.write(fn * ".csv", df2; decimal = ',', delim='\t')
        XLSX.writetable(fn * ".xlsx", collect(DataFrames.eachcol(df2)), DataFrames.names(df2), overwrite = true)

    end

end
