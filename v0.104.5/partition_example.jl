
using Oceananigans
using Oceananigans.DistributedComputations: Equal
using MPI
MPI.Init()

partition = Partition(y=Equal())

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @show partition
end
