@doc raw"""
	long_duration_storage_sparse_chronology!(EP::Model, inputs::Dict, setup::Dict)
This function creates variables and constraints enabling modeling of long duration storage resources when modeling representative time periods using the [sparse chronology strategy](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5061243).\
"""
function long_duration_storage_sparse_chronology!(EP::Model, inputs::Dict, setup::Dict)
    println("Long Duration Storage Sparse Chronology Module")

    gen = inputs["RESOURCES"]

    CapacityReserveMargin = setup["CapacityReserveMargin"]

    STOR_LONG_DURATION = inputs["STOR_LONG_DURATION_SPARSE_CHRONOLOGY"]
    hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod
    dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
    NPeriods = size(dfPeriodMap)[1] # Number of modeled periods
    println("NPeriods: ",NPeriods)

    REP_PERIODS = unique(dfPeriodMap[!,"Rep_Period"])
    NRepPeriods = size(REP_PERIODS)[1]
    println("NRepPeriods: ",NRepPeriods)
    println("REP_PERIODS      : ",REP_PERIODS)

    # compute the length of each "partition" where a partition is a consequtive set of repeated representative periods
    PARTITION_ENDS_MASK = [dfPeriodMap[1:end-1,"Rep_Period"] .!= dfPeriodMap[2:end,"Rep_Period"]; true]
    PARTITION_ENDS = findall(PARTITION_ENDS_MASK)
    PARTITION_LENGTHS = PARTITION_ENDS - [0; PARTITION_ENDS[1:end-1]]
    PARTITION_REP_PERIOD_INDEX = dfPeriodMap[PARTITION_ENDS_MASK,"Rep_Period_Index"]
    println("PARTITION_ENDS = ",PARTITION_ENDS)
    println("PARTITION_LENGTHS = ",PARTITION_LENGTHS)
    println("sum(PARTITION_LENGTHS) = ",sum(PARTITION_LENGTHS))
    println("PARTITION_REP_PERIOD_INDEX = ",PARTITION_REP_PERIOD_INDEX)
    NPartitions = size(PARTITION_LENGTHS)[1]
    println("NPartitions: ",NPartitions)
    PARTITIONS_LONGER_THAN_ONE = [p for p in 1:NPartitions if PARTITION_LENGTHS[p]>1]
    println("PARTITIONS_LONGER_THAN_ONE = ",PARTITIONS_LONGER_THAN_ONE)
    REP_PERIOD_INDEX_PARTITIONS_LONGER_THAN_ONE = sort(unique(PARTITION_REP_PERIOD_INDEX[PARTITIONS_LONGER_THAN_ONE]))
    println("REP_PERIOD_INDEX_PARTITIONS_LONGER_THAN_ONE = ",REP_PERIOD_INDEX_PARTITIONS_LONGER_THAN_ONE)

    eTotalCapEnergy = EP[:eTotalCapEnergy]
    vCHARGE = EP[:vCHARGE] # charge power
    vS = EP[:vS] # stored energy
    vP = EP[:vP] # discharge power

    ### Implement sparse chonology equations and constraints from https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5061243 ###
    # Total Constraints = (NRepPeriods + 3*NPartitions + 2*NPartitionsLongerThanOne + 2*NRepPeriods*(hours_per_subperiod-1))*NBatteries
    
    # (10): change in energy store over representative period
    @variable(EP,   vDeltaStoreInRepPeriod[y in STOR_LONG_DURATION, r=1:NRepPeriods])
    @constraint(EP, cDeltaStoreInRepPeriod[r=1:NRepPeriods, y in STOR_LONG_DURATION], # cSoCBalLongDurationStorageStart. Need to index by r,y rather than y,r due to array division in write_storagedual.jl
                vDeltaStoreInRepPeriod[y,r] == vS[y, hours_per_subperiod*r] - vS[y, hours_per_subperiod*(r-1)+1]) # NRepPeriods*NBatteries

    # (11): Connect initial SOC for each partition to the next partition
    @variable(EP,       vPartitionInitialStore[y in STOR_LONG_DURATION, p=1:NPartitions] >= 0)
    @constraint(EP, cLinkPartitionInitialStore[y in STOR_LONG_DURATION, p=1:NPartitions], # NPartitions*NBatteries
                vPartitionInitialStore[y, mod1(p+1, NPartitions)] == vPartitionInitialStore[y,p] + PARTITION_LENGTHS[p]*vDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]] ) # self-discharge term?
    
    # (12): Define a variable which tracks the maximum increase in SOC within each representative period
    @variable(EP,   vMaxDeltaStoreInRepPeriod[y in STOR_LONG_DURATION, r=1:NRepPeriods] >= 0)
    @constraint(EP, cMaxDeltaStoreInRepPeriod[y in STOR_LONG_DURATION, r=1:NRepPeriods, h=2:hours_per_subperiod], # NRepPeriods*(hours_per_subperiod-1)*NBatteries
                vMaxDeltaStoreInRepPeriod[y,r] >= vS[y, hours_per_subperiod*(r-1)+h] - vS[y, hours_per_subperiod*(r-1)+1] )

    # (13): Define a variable which tracks the maximum decrease (minimum over changes) in SOC within each representative period
    @variable(EP,   vMinDeltaStoreInRepPeriod[y in STOR_LONG_DURATION, r=1:NRepPeriods] <= 0)
    @constraint(EP, cMinDeltaStoreInRepPeriod[y in STOR_LONG_DURATION, r=1:NRepPeriods, h=2:hours_per_subperiod], # NRepPeriods*(hours_per_subperiod-1)*NBatteries
                vMinDeltaStoreInRepPeriod[y,r] <= vS[y, hours_per_subperiod*(r-1)+h] - vS[y, hours_per_subperiod*(r-1)+1] )

    # (14): Constrain minimum SOC in first rep period in each partition to be greater than zero
    @constraint(EP, cPartitionMinStoreInFirstRepPeriod[y in STOR_LONG_DURATION, p=1:NPartitions], # NPartitions*NBatteries
                0                  <= vPartitionInitialStore[y,p] + vMinDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]])

    # (15): Constrain minimum SOC in last rep period in each partition to be greater than zero (only needed if this partition is longer than one period!)
    @constraint(EP, cPartitionMinStoreInLastRepPeriod[y in STOR_LONG_DURATION, p=PARTITIONS_LONGER_THAN_ONE], # NPartitionsLongerThanOne*NBatteries
                0                  <= vPartitionInitialStore[y,p] + vMinDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]] + (PARTITION_LENGTHS[p]-1)*vDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]])

    # (16): Constrain maximum SOC in first rep period in each partition to be less than the total capacity
    @constraint(EP, cPartitionMaxStoreInFirstRepPeriod[y in STOR_LONG_DURATION, p=1:NPartitions], # NPartitions*NBatteries
                eTotalCapEnergy[y] >= vPartitionInitialStore[y,p] + vMaxDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]])

    # (17): Constrain maximum SOC in last rep period in each partition to be less than the total capacity (only needed if this partition is longer than one period!)
    @constraint(EP, cPartitionMaxStoreInLastRepPeriod[y in STOR_LONG_DURATION, p=PARTITIONS_LONGER_THAN_ONE], # NPartitionsLongerThanOne*NBatteries
                eTotalCapEnergy[y] >= vPartitionInitialStore[y,p] + vMaxDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]] + (PARTITION_LENGTHS[p]-1)*vDeltaStoreInRepPeriod[y, PARTITION_REP_PERIOD_INDEX[p]])


    # vS at start of each representative period is a free parameter in principle. In practice it effects capacity reserve constraints in storage_all.jl
    @constraint(EP, cSoCBalLongDurationStorageStart[y in STOR_LONG_DURATION, r=1:NRepPeriods],
                vS[y, hours_per_subperiod*(r-1)+1] == eTotalCapEnergy[y]/2)

end

