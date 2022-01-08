module ContractMachines

#-- Required modules
# Wiring diagrams
using AlgebraicDynamics
using AlgebraicDynamics.DWDDynam
using Catlab.WiringDiagrams
import AlgebraicDynamics: oapply

# ode solver
using DifferentialEquations
import DifferentialEquations: ODEProblem

# Contract display
using IntervalSets
using PrettyTables
using Printf

# Static contracts
include("StaticContracts.jl")
using .StaticContracts	# Must have period as submodule

#-- Information accesible by using module
export ContractMachine, oapply, check_contract, ODEProblem

#-- Datatypes used by contract machine
const ContractOutputTable = Dict{KEY, NamedTuple{ (:input, :output),
                                                  Tuple{ Vector{ Pair{Symbol, Bool} },
                                                         Vector{ Pair{Symbol, Bool} } } }
                                 } where KEY <: Union{Symbol, Any}

const ContractTimeTable = Dict{Any, Any}

const ContractOutputBox = NamedTuple{ (:input, :output),
                                      Tuple{ Vector{Bool}, Vector{Bool} } }

#-- Type to display pretty table
struct ContractTable
    table::Union{ContractOutputTable, ContractTimeTable, ContractOutputBox}
end

#-- Contracts are defined via intervals
struct ContractMachine{T<:Real}
    static_contract::StaticContract{T}
    machine::AbstractMachine{T}
    fcontract::Function

    # inner constructor
    function ContractMachine{T}(static_contract::StaticContract{T}, machine::AbstractMachine{T};
                                fcontract = nothing) where T<:Real
        # contract function
        # Only define function in none is provided, used to not overwrite the composed function from oapply
        if fcontract == nothing
            fcontract = (u::AbstractVector, x::AbstractVector, p=nothing, t=0) -> begin
                # check whether contracts at input ports are satisfied
                Rin = map( (xin,cont) -> xin in cont, x, static_contract.input )

                # check whether contracts at output ports are satisfied
                Rout = map( (rout,cont) -> rout in cont, readout(machine)(u, p, t), static_contract.output )
                return ContractTable( (input = Rin, output = Rout) )
            end
        end
        # make a new machine satisfying these restrictions
        new{T}(static_contract, machine, fcontract)
    end
end

# outer constructor
function ContractMachine{T}(cinput::Vector, nstates::Int, coutput::Vector,
                            dynamics::Function, readout::Function;
                            fcontract = nothing, mtype = :continuous) where T<:Real
    # contracts
    static_contract = StaticContract{T}(cinput, coutput)

    # each element in the vectors is a contract for the port on a box
    ninputs = length(cinput)
    noutputs = length(coutput)

    # select a machine type
    if mtype == :continuous
        machine = ContinuousMachine{T}(ninputs, nstates, noutputs, dynamics, readout)
    elseif mtype :discrete
        machine = DiscreteMachine{T}(ninputs, nstates, noutputs, dynamics, readout)
    end
    # make a machine using inner constructor
    ContractMachine{T}(static_contract, machine, fcontract=fcontract )
end

#-- Compose multiple contract machines --

function oapply(d::WiringDiagram, ms::Vector{ContractMachine{T}}) where T<:Real
    # compose contracts
    static_contract = StaticContracts.oapply(d, map(m -> m.static_contract, ms))

    # Get the initial index of the state vector of each each box.
    # The composition concatenates these vectors into a single column.
    nstate = 0
    index = Array{ UnitRange{Int} }(undef, nboxes(d))

    for i in 1:nboxes(d)
        nstate += nstates(ms[i].machine)                           # initial index of vector
        index[i] = nstate - nstates(ms[i].machine) + 1 : nstate    # store vector indeces
    end

    # get the name of each box
    name = map( box -> box.value, boxes(d) )

    #---- Evaluate contracts
    function fcontract(u::AbstractVector, x::AbstractVector, p=nothing, t=0)
        # get the output of each box
        rout = map( id -> readout(ms[id].machine)( u[index[id]], p, t ), 1:nboxes(d) )

        # evaluate the contract function
        fout = Array{Dict}(undef, nboxes(d))

        for id in 1:nboxes(d)                       # check all boxes in the diagram
            # collect the inputs of each box
            xin = zeros( ninputs(ms[id].machine) )

            for w in in_wires(d, id)                # check all wires going into a box
                if w.source.box != input_id(d)      # iternal inputs use are the readout of another box
                    xin[w.target.port] += rout[w.source.box][w.source.port]
                else                                # external inputs make use of a given vector
                    xin[w.target.port] += x[w.source.port]
                end
            end

            # evaluate the contract function of each box
            foutput = ms[id].fcontract( u[index[id]], xin, p, t ).table

            # for atomic boxes, assign the box name to the inputs and outputs
            if typeof(foutput) <: NamedTuple
                fout[id] = Dict(name[id] => (input = input_ports(d, id) .=> foutput.input,
                                             output = output_ports(d, id) .=> foutput.output))
            else    # for nested boxes, append the name of the parent box to the child boxes
                fout[id] = Dict( (name[id] .=> keys(foutput)) .=> values(foutput) )
            end
        end
        # return a contract table with a box directory assigned to inputs and outputs
        return ContractTable( merge(fout...) )
    end

    # compose the dynamics functions
    machine = DWDDynam.oapply(d, map(m -> m.machine, ms))

    # return the composed machine
    return ContractMachine{T}(static_contract, machine, fcontract=fcontract)
end

# Identify during which time intervals a signal is zero [false == contract is violated]
function failureInterval(arr::AbstractVector, time=nothing, out_type="time")
    # check each element in the array
    index = Array{Tuple}(undef, 0)
    start = 1

    for i in 1:length(arr)             # no previous element
        if i != 1
            if arr[i] > arr[i-1]       # rising signal
                push!( index, (start, i-1) )

            elseif arr[i] < arr[i-1]   # falling signal
                start = i
            end

            if i == length(arr) && arr[end] <= 0   # end of interval
                push!( index, (start, i) )
            end
        end
    end

    # map indeces to given array
    if time != nothing && out_type == "time"
        index = map( set -> (time[set[1]], time[set[2]]), index)
    end

    # output 2d array of start and stop times
    return index
end

# Assign contract failure times to each wire of each box
function check_contract(sol::T1, machine::ContractMachine{T2}, x0::AbstractVector, p=nothing, t=0;
                        out_type="time" ) where {T1<:ODESolution, T2<:Real}

    # evalutate the contract function throughout time interval
    fout = map( t -> machine.fcontract(sol(t), x0, p, t), sol.t )

    # Store time at which each contract fails for each box directory:
    # Initial evaluation of contracts. Used to get the names and wires of each box.
    set0 = fout[begin].table
    dict = Dict()

    # go through all directories
    for key in keys(set0)
        # function to map indeces of contract outputs
        mapOut = index -> map(i -> set0[key][index][i].first => begin
                              # array of the contract state for each time step
                              contract = map(set -> set.table[key][index][i].second, fout);
                              # find the duration of 0's in the array
                              failureInterval(contract, sol.t, out_type);
                              end,
                              1:length(set0[key][index]) )
        # store the times at the directory
        dict[key] = ( input=mapOut(:input), output=mapOut(:output) )
    end

    # output a data structure like that of fcontract but with failure intervals.
    return ContractTable(dict)
end

#-- Pretty printing of output

# display contract machines as product of intervals
function Base.show(io::IO, vf::ContractMachine)
    display(vf.static_contract)
end

# display contract tables (output of fcontract and eval_contract) as a pretty table
function Base.show(io::IO, data::ContractTable)
    # get the data structure
    dict = data.table

    # select the appriate table
    if typeof(dict) <: ContractOutputBox        # atomic box
        # assign a contract state to each port
        mapOut = index -> join( string.(1:length(dict[index])) .* " : " .* string.(dict[index]), "\n" )

        # do not assign a box as there is no diagram
        output = hcat( [mapOut(:input)], [mapOut(:output)] )
        header = ( ["input", "output"], ["port: contract", "port: contract"] )

    else # composed boxes
        box = string.(keys(dict))

        if typeof(dict) <: ContractOutputTable  # output of fcontract
            mapOut = index -> map( x -> join(
                map(pair -> string(pair.first) * " : " * string(pair.second), x[index]), "\n"),
                                   values(dict) )

            header =  ( ["box", "input", "output"], ["directory", "wire: contract", "wire: contract"] )

        else    # output of eval_contract
            mapOut = index::Symbol -> map( value -> begin
                                           # make array with strings of time intervals for each box
                                           # Note: need way to change formatting for intergers
                                           out = map( pair ->
                                                      map( set -> @sprintf("%s : %f , %f", pair.first, set...), pair.second ),
                                                      value[index] );

                                           # combine strings into large string seperated by line jumps
                                           join( vcat(out...), "\n")
                                           end,
                                           values(dict) )

            header = ( ["box", "input", "output"], ["directory", "wire: failure interval", "wire: failure interval"] )
        end
        # both types have the same structure
        output = hcat(box, mapOut(:input), mapOut(:output))
    end
    # display a pretty table
    pretty_table(output, header = header, backend=:html, linebreaks=true, tf=tf_html_minimalist, alignment=:c)
end

#-- helper functions

# compose machines given the name of each box
function oapply(d::WiringDiagram, ms::Dict{Symbol, ContractMachine{T}}) where T<:Real
    oapply(d, map(box -> ms[box.value], boxes(d)) )
end

# solve the dynamics of the contract machines
function ODEProblem(m::ContractMachine{T}, u0::AbstractVector, xs::AbstractVector, tspan::Tuple, p=nothing) where T<:Real
    ODEProblem(m.machine, u0, xs, tspan, p)
end

end # module
