function write_opwrap_lds_stor_init(path::AbstractString,
        inputs::Dict,
        setup::Dict,
        EP::Model)
    ## Extract data frames from input dictionary
    gen = inputs["RESOURCES"]
    zones = zone_id.(gen)

    G = inputs["G"]

    # Initial level of storage in each modeled period
    NPeriods = size(inputs["Period_Map"])[1]
    dfStorageInit = DataFrame(Resource = inputs["RESOURCE_NAMES"], Zone = zones)
    socw = zeros(G, NPeriods)
    for i in 1:G
        if i in inputs["STOR_LONG_DURATION"]
            socw[i, :] = value.(EP[:vSOCw])[i, :]
        end
        # if i in inputs["STOR_LONG_DURATION_SPARSE_CHRONOLOGY"]
        #     socw[i, :] = value.(EP[:vSOCw])[i, :]
        # end
        if !isempty(inputs["VRE_STOR"])
            if i in inputs["VS_LDS"]
                socw[i, :] = value.(EP[:vSOCw_VRE_STOR][i, :])
            end
        end
    end
    if setup["ParameterScale"] == 1
        socw *= ModelScalingFactor
    end

    dfStorageInit = hcat(dfStorageInit, DataFrame(socw, :auto))
    auxNew_Names = [Symbol("Resource"); Symbol("Zone"); [Symbol("n$t") for t in 1:NPeriods]]
    rename!(dfStorageInit, auxNew_Names)
    CSV.write(joinpath(path, "StorageInit.csv"),
        dftranspose(dfStorageInit, false),
        header = false)

    # Write storage evolution over full time horizon
    hours_per_subperiod = inputs["hours_per_subperiod"];
    t_interior = 2:hours_per_subperiod
    T_hor = hours_per_subperiod*NPeriods # total number of time steps in time horizon
    SOC_t = zeros(G, T_hor)
    stor_lds = inputs["STOR_LONG_DURATION"]
    stor_lds_sc = inputs["STOR_LONG_DURATION_SPARSE_CHRONOLOGY"]
    stor_hydro_lds = inputs["STOR_HYDRO_LONG_DURATION"]
    period_map = inputs["Period_Map"].Rep_Period_Index
    pP_max = inputs["pP_Max"]
    e_total_cap = value.(EP[:eTotalCap])
    v_charge = value.(EP[:vCHARGE])
    v_P = value.(EP[:vP])
    if setup["ParameterScale"] == 1
        v_charge *= ModelScalingFactor
        v_P *= ModelScalingFactor
    end
    if !isempty(stor_hydro_lds)
        v_spill = value.(EP[:vSPILL])
    end
    for r in 1:NPeriods
        w = period_map[r]
        t_r = hours_per_subperiod * (r - 1) + 1
        t_start_w = hours_per_subperiod * (w - 1) + 1
        t_interior = 2:hours_per_subperiod

        if !isempty(stor_lds)
            SOC_t[stor_lds, t_r] = socw[stor_lds, r] .* (1 .- self_discharge.(gen[stor_lds])) .+ efficiency_up.(gen[stor_lds]) .* v_charge[stor_lds, t_start_w] .- 1 ./ efficiency_down.(gen[stor_lds]) .* v_P[stor_lds, t_start_w]

            for t_int in t_interior
                t = hours_per_subperiod * (w - 1) + t_int
                SOC_t[stor_lds, t_r + t_int - 1] = SOC_t[stor_lds, t_r + t_int - 2] .* (1 .- self_discharge.(gen[stor_lds])) .+ efficiency_up.(gen[stor_lds]) .* v_charge[stor_lds, t] .- 1 ./ efficiency_down.(gen[stor_lds]) .* v_P[stor_lds, t]
            end
        end

        if !isempty(stor_hydro_lds)
            SOC_t[stor_hydro_lds, t_r] = socw[stor_hydro_lds, r] .- 1 ./ efficiency_down.(gen[stor_lds]) .* v_P[stor_hydro_lds, t_start_w] .- v_spill[stor_hydro_lds, t_start_w] .+ pP_max[stor_hydro_lds, t_start_w] .* e_total_cap[stor_hydro_lds]

            for t_int in t_interior
                t = hours_per_subperiod * (w - 1) + t_int
                SOC_t[stor_hydro_lds, t_r + t_int - 1] = SOC_t[stor_hydro_lds, t_r + t_int - 2] .- 1 ./ efficiency_down.(gen[stor_hydro_lds]) .* v_P[stor_hydro_lds, t] .- v_spill[stor_hydro_lds, t] .+ pP_max[stor_hydro_lds, t] .* e_total_cap[stor_hydro_lds]
            end
        end
    end

    if !isempty(stor_lds_sc)
        dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
        # compute the length of each "partition" where a partition is a consequtive set of repeated representative periods
        PARTITION_ENDS_MASK = [dfPeriodMap[1:end-1,"Rep_Period"] .!= dfPeriodMap[2:end,"Rep_Period"]; true]
        PARTITION_ENDS = findall(PARTITION_ENDS_MASK)
        PARTITION_STARTS = [0; PARTITION_ENDS[1:end-1]] .+ 1
        PARTITION_LENGTHS = PARTITION_ENDS - PARTITION_STARTS .+ 1
        PARTITION_REP_PERIOD_INDEX = dfPeriodMap[PARTITION_ENDS_MASK,"Rep_Period_Index"]
        NPartitions = size(PARTITION_LENGTHS)[1]

        PARTITION_STARTS_HOUR = 1 .+ (PARTITION_STARTS.-1)*hours_per_subperiod

        vPartitionInitialStore = value.(EP[:vPartitionInitialStore])[stor_lds_sc, :]
        if setup["ParameterScale"] == 1
            vPartitionInitialStore *= ModelScalingFactor
        end
        SOC_t[stor_lds_sc, PARTITION_STARTS_HOUR] = vPartitionInitialStore

        PARTITION_LENGTHS_HOUR = PARTITION_LENGTHS*hours_per_subperiod
        PARTITION_REP_PERIOD_INDEX_HOUR = 1 .+ (PARTITION_REP_PERIOD_INDEX.-1)*hours_per_subperiod
        for p in 1:NPartitions
            t_interior   = 1:PARTITION_LENGTHS_HOUR[p]-1
            t_partition  = PARTITION_STARTS_HOUR[p].+t_interior.-1
            t_rep_period = PARTITION_REP_PERIOD_INDEX_HOUR[p]:PARTITION_REP_PERIOD_INDEX_HOUR[p].+hours_per_subperiod.-1

            for t_i in t_interior
                t_p = t_partition[t_i]
                t_r = t_rep_period[mod1(t_i, hours_per_subperiod)]
                SOC_t[stor_lds_sc, t_p+1] = SOC_t[stor_lds_sc, t_p] .+ efficiency_up.(gen[stor_lds_sc]).*v_charge[stor_lds_sc, t_r] .- 1 ./ efficiency_down.(gen[stor_lds_sc]).*v_P[stor_lds_sc, t_r]
            end
        end
    end
    
    df_SOC_t = DataFrame(Resource = inputs["RESOURCE_NAMES"], Zone = zones)
    df_SOC_t = hcat(df_SOC_t, DataFrame(SOC_t, :auto))
    auxNew_Names = [Symbol("Resource"); Symbol("Zone"); [Symbol("n$t") for t in 1:T_hor]]
    rename!(df_SOC_t,auxNew_Names)
    CSV.write(joinpath(path, "StorageEvol.csv"), dftranspose(df_SOC_t, false), writeheader=false)

end
