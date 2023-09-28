

function parse_numbers(path::AbstractString)
    numbers = Vector{Int32}()
    for node in split(path, " ")
        number = parse(Int32, view(node, 1:length(node)-1))
        # number = node[end] == "-" ? flipnode(number) : noflip(number)
        push!(numbers, number)
    end 
    return numbers
end

function get_stats_gfa2(f::String)
    h = open(f, "r")
    unique_nodes = Set{Int}()
    edges = 0
    for line in eachline(h)
        identifier, path = split(line, "\t")
        node_ids = parse_numbers(path)
        for node in node_ids
            push!(unique_nodes, node)
        end
        edges += length(edges)-1 
    end
    println("Nodes: ", nodes)
    println("Edges: ", edges)
end

get_stats_gfa2("/net/phage/linuxhome/mgx/people/rickb/Campy/complete_genomees/graph/complete_campy_graph_gfa.gfa2")