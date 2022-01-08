module StaticContracts

#-- Modules
# wiring diagrams
using Catlab.WiringDiagrams
import Catlab.WiringDiagrams: oapply  # Needed to overwrite oapply function

# intervals of real numbers
using IntervalSets

#-- Accessible functions
export StaticContract, oapply

#-- Datatype that accepts all bound conditions
const NInterval = Vector{ Union{ClosedInterval{T},
                                OpenInterval{T},
                                Interval{:open, :closed, T},
                                Interval{:closed, :open, T} } } where T<:Real

#-- Contracts are defined via intervals:
struct StaticContract{T<:Real}
    ninputs::Int
    noutputs::Int
    input::NInterval{T}
    output::NInterval{T}

    # Constructor
    function StaticContract{T}(input::Vector, output::Vector) where T<:Real
        ninputs = length(input)
        noutputs = length(output)

        # cannot have empty contract
        for contract in [input; output]
            if isempty(contract)
                error("the interval $contract is backwards")
            end
        end

        new{T}(ninputs, noutputs, input, output)
    end
end

#-- Compose contracts with diagram: DWDDynam.
function oapply(d::WiringDiagram, ms::Vector{StaticContract{T}}) where T
    # boxes in diagram
    box = boxes(d)

    # Ensure machines fill diagram
    if nboxes(d) != length(ms)
        error("there are $nboxes(d) boxes but $length(ms) machines")
    end
    # each wire must have a contract
    for id in 1:nboxes(d)
        b = box[id]
        if length(b.input_ports) != ms[id].ninputs || length(b.output_ports) != ms[id].noutputs
            name = b.value
            error("number of ports do not match number of contracts at box $name")
        end
    end

    # Check all contracts inside the diagram (excluding wires that enter or exit the diagram)
    for w in wires(d, :Wire)
        # Variable names, see section 3.3.3, P.16
        target_var = box[w.target.box].input_ports[w.target.port]
        source_var = box[w.source.box].output_ports[w.source.port]

        # ensure target and source name match.
        if target_var != source_var
            error("variable names do not match at $w")
        end

        # Induced contract, see section 3.3.1, P.13
        # box names
        source_name = box[w.source.box].value
        target_name = box[w.target.box].value

        # check if source contract is compatible with target contract
        contract_source = ms[w.source.box].output[w.source.port]
        contract_target = ms[w.target.box].input[w.target.port]
        overlap = intersect(contract_source, contract_target)

        if isempty(overlap) == true
            error("contract $contract_source of $source_name does not satisfy contract $contract_target of $target_name at wire $target_var")
        elseif overlap != contract_source
            println("WARNING: contract $contract_source of $source_name lies outside contract $contract_target of $target_name at wire $target_var")
        end
    end

    # Assign external contracts
    input = map( w -> ms[w.target.box].input[w.target.port], wires(d, :InWire))
    output = map( w -> ms[w.source.box].output[w.source.port], wires(d, :OutWire))

    # Create new machine
    return StaticContract{T}(input, output)
end

# Compose using a library DWDDynam.
function oapply(d::WiringDiagram, ms::Dict{Symbol, StaticContract{T}}) where T
    oapply(d, map(box -> ms[box.value], boxes(d)) )
end

#-- Display struct as product of intervals:
function Base.show(io::IO, vf::StaticContract)
    # Check all contracts
    list = [vf.input; vf.output]
    output = ""

    # iterate across all contracts
    len = length(list)

    for i in 1:len
        contract = list[i]
        # Switch symbol for infinity
        if -contract.left == contract.right == Inf
            output *= "ℝ"
        else
            left = contract.left == -Inf ? "-∞" : contract.left
            right = contract.right == Inf ? "∞"  : contract.right
            output *= "[$left,$right]"
        end
        # add product operator [spaces make it easier to read]
        if i != len
            output *= " × "
        end
    end

    # display combined string
    print("StaticContract( $output )")
end

end # module
