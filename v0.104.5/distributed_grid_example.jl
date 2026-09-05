
using Oceananigans
using Oceananigans.DistributedComputations
using MPI; MPI.Init()

child_architecture = CPU()
architecture = Distributed(child_architecture)

grid = RectilinearGrid(architecture,
                       size = (48, 48, 16),
                       x = (0, 64),
                       y = (0, 64),
                       z = (0, 16),
                       topology = (Periodic, Periodic, Bounded))

@handshake @info grid
