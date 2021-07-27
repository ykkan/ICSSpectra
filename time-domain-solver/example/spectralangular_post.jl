using MPI
using DelimitedFiles
using FFTW
########################
#### user parameters setting
#######################

const fileprefix = "radiation"
const charge = -1.0
const lambda0 = 1000e-9

###########################
const dirpath = isempty(ARGS) ? "./" : ARGS[1]
######################
### parse input file
detectorinfo_raw = joinpath(dirpath, "detectorinfo.dat") |> readlines
thxrange_raw = split(detectorinfo_raw[1], " ")
thyrange_raw = split(detectorinfo_raw[2], " ")
pulserange_raw = split(detectorinfo_raw[3], " ")

const thx_min = parse(Float64, thxrange_raw[2])
const thx_max = parse(Float64, thxrange_raw[3])
const nthx = parse(Int, thxrange_raw[4])

const thy_min = parse(Float64, thyrange_raw[2])
const thy_max = parse(Float64, thyrange_raw[3])
const nthy = parse(Int, thyrange_raw[4])


const tu_min = parse(Float64, pulserange_raw[2])
const tu_max = parse(Float64, pulserange_raw[3])
const ntsample = parse(Int, pulserange_raw[4])


#### physics constant
c0 = 3.0e8
q0 = 1.602e-19
eps0 = 8.854e-12
hbar = 1.054571e-34

### derived parameters
const nth = nthx != 1 ? nthx : nthy
const dth = nthx != 1 ?  ( thx_max - thx_min ) / (nthx - 1) : ( thy_max - thy_min ) / (nthy - 1)
const (th_min, th_max) = nthx != 1 ? (thx_min, thx_max) : (thy_min, thy_max)
const Tu = tu_max - tu_min

const w0 = 2 * pi * c0 / lambda0
const dtu = Tu / ntsample

const nwsample = ceil(Int, ntsample / 2.0) # inlcuding zero frequency
const amp = charge^2 * q0^2 * dtu^2 / (16.0 * pi^3 * eps0 * c0)



### MPI initialization
MPI.Init()
comm = MPI.COMM_WORLD
nprocs = MPI.Comm_size(comm)
myid = MPI.Comm_rank(comm)
rootid = 0
######### main routine
nth_local_list = fill(floor(Int, nth / nprocs), nprocs)
nth_local_list[1:nth % nprocs] .+= 1
nth_local = nth_local_list[myid + 1]

if myid == 0
    startindex_global = 1
else
    startindex_global = 1 + sum(nth_local_list[1:myid])
end
endindex_global = startindex_global + nth_local_list[myid + 1] - 1



angularspectrums_local = zeros(nwsample, nth_local)
begin
    for i = 1:nth_local
        filepath = "rawdata/" * fileprefix * "_" * string(startindex_global + i - 1) * ".dat"
        angularspectrums_local[:, i] = amp * sum(abs2.(mapslices(fft, readdlm(filepath), dims=1)), dims=2)[1:nwsample]
    end
end

dwbar = 2 * pi / Tu;
wbars = dwbar * [0:(nwsample - 1);]
ws = w0 * wbars
dw = w0 * dwbar
th_min_local = -dth * (startindex_global - ceil(nth / 2) )
th_max_local =  dth * (  endindex_global - ceil(nth / 2) )

nwsample_dump = nwsample
angularspectrums_dump = zeros(nwsample_dump, nth)
angularspectrums_dump[:, startindex_global:endindex_global] = angularspectrums_local[1:nwsample_dump, :]
angularspectrums_dump = MPI.Reduce(angularspectrums_dump, +, rootid, comm)
if myid == rootid
    open("spectralangular.dat", "w") do io
        writedlm(io, angularspectrums_dump)
    end
end
