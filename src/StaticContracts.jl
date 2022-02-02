module StaticContracts

#-- Modules
# wiring diagrams
using Catlab.WiringDiagrams
import Catlab.WiringDiagrams: oapply  # Needed to overwrite oapply function

# intervals of real numbers
using IntervalSets

#-- Accessible functions
export StaticContract, oapply, StaticContractDict

#-- Contracts are defined via intervals:
struct StaticContract{T<:Real}
    input::Vector{Interval}
    output::Vector{Interval}

    # Constructor
    function StaticContract{T}(input::Vector, output::Vector) where T<:Real
        # cannot have empty contract
        for contract in [input; output]
            if isempty(contract)
                error("the interval $contract is empty or backwards")
            end
        end
        new{T}(input, output)
    end
end


#-- Display struct as product of intervals:
function formatInterval(contract::Interval) 
    # Switch symbol for infinity
    if -contract.left == contract.right == Inf
        output = "ℝ"
    else
        # check left bound
        if( contract.left == -Inf )
            left = "-∞"
            lbracket = "("
        else
            left = contract.left
            lbracket = isleftopen(contract) ? "(" : "["
        end

        # check right bound
        if( contract.right == Inf )
            right = "∞"
            rbracket = ")"
        else
            right = contract.right
            rbracket = isrightopen(contract) ? ")" : "]"
        end    

        # Concatenate strings
        output = "$lbracket" * "$left,$right" * "$rbracket"
    end
    return output
end

function Base.show(io::IO, vf::StaticContract)
    # Check all contracts
    list = [vf.input; vf.output]
    output = ""

    # iterate across all contracts
    len = length(list)

    for i in 1:len
        contract = list[i]
        output *= formatInterval(contract)
        # add product operator [spaces make it easier to read]
        if i != len
            output *= " × "
        end
    end
    # display combined string
    print("StaticContract( $output )")
end

#-- Compose contracts with diagram: 
function oapply(d::WiringDiagram, ms::Vector{StaticContract{T}}) where T
    # boxes in diagram
    box = boxes(d)

    # Ensure machines fill diagram
    n_contract = length(ms)
    
    if nboxes(d) != n_contract
        error("there are $nboxes(d) boxes but $length(ms) machines")
    end
    
    # Boxes must have unique names to be identified
    if n_contract > 1 
        name = map( m -> m.value, box )
        for i in 1:(n_contract-1)
            curr_name = name[i]
            if curr_name in name[ (i+1):n_contract ]
                error("Two or more boxes share the name: $curr_name")
            end
        end
    end
    
    # Each wire must have a contract 
    for id in 1:nboxes(d)
        curr_box = box[id]
        if length(curr_box.input_ports) != length(ms[id].input) || length(curr_box.output_ports) != length(ms[id].output)
            name = curr_box.value
            error("number of ports do not match number of contracts at $name (id=$id)")
        end
    end

    # Check all contracts inside the diagram (excluding wires that enter or exit the diagram)
    for w in wires(d, :Wire)
        # id number of boxes
        source_id = w.source.box
        target_id = w.target.box
        # box names
        source_name = box[source_id].value
        target_name = box[target_id].value
        
        # 1. Variable names, see section 3.3.3, P.16
        source_var = box[source_id].output_ports[w.source.port]
        target_var = box[target_id].input_ports[w.target.port]
        
            # ensure target and source name match.
        if target_var != source_var
            error("wire \"$target_var\" of $source_name (id=$source_id) " * 
                  "does not match wire \"$source_var\" of $target_name (id=$target_id)")
        end
        
        # 2. Induced contract, see section 3.3.1, P.13
        contract_source = ms[source_id].output[w.source.port] 
        contract_target = ms[target_id].input[w.target.port]
        overlap = intersect(contract_source, contract_target)
        
        format_source = formatInterval(contract_source)
        format_target = formatInterval(contract_target)

            # check if source contract is compatible with target contract
        if isempty(overlap) == true
            error("Incompatible contract between $source_name (id=$source_id) " * 
                  "and $target_name (id=$target_id) at wire \"$target_var\": " * 
                  "$format_source ∩ $format_target = ∅")
            
            # Check for undefined behavior
        elseif overlap != contract_source
            println("Intervals do not overlap between $source_name (id=$source_id) " * 
                    "and $target_name (id=$target_id) at wire \"$target_var\": " *
                    "$format_source ∩ $format_target ≠ $format_source")
        end
    end

    # Assign external contracts
    input = map( w -> ms[w.target.box].input[w.target.port], wires(d, :InWire))
    output = map( w -> ms[w.source.box].output[w.source.port], wires(d, :OutWire))

    # Create new machine
    return StaticContract{T}(input, output)
end

# Compose using a library 
const StaticContractDict{T} = Dict{ Symbol, StaticContract{T} }

function oapply(d::WiringDiagram, ms::StaticContractDict{T}) where T
    box = boxes(d)
    name_contract = keys(ms)
    
    # add unassigned contracts
    if length(box) != length(name_contract)
        for b in box
            name = b.value
            if !(name in name_contract)
                # create default box
                n_input = length(b.input_ports)
                n_output = length(b.output_ports)
                default_contract = StaticContract{T}(repeat([-Inf..Inf], n_input), 
                                                     repeat([-Inf..Inf], n_output) )
                # add to dictionary
                ms[name] = default_contract
            end
        end
    end
    # compose with array 
    oapply(d, map(box -> ms[box.value], box) )
end

end # module
