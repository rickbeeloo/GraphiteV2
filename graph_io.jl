function parse_numbers(path::AbstractString)
    numbers = Vector{Int32}()
    for node in split(path, " ")
        number = parse(Int32, view(node, 1:length(node)-1))
        number = node[end] == "-" ? flipnode(number) : noflip(number)
        push!(numbers, number)
    end 
    return numbers
end

function read_queries(f::String, query_ids::OrderedSet{String})
    queries = Vector{Vector{Int32}}()
    for (i, line) in enumerate(eachline(f))
        identifier, path = split(line, "\t")
        if identifier in query_ids
            path_numbers = parse_numbers(path)
            push!(query_ids, identifier)
            push!(queries, path_numbers)
        end
    end
    return queries
end

function read_node_sizes(f::String)
    size_map = Dict{Int32, Int32}()
    for line in eachline(f)
        node_id, seq = split(line)
        seq_size = Int32(length(seq))
        node_id = parse(Int32, node_id)
        size_map[node_id] = seq_size
    end
    return size_map
end

function processGFA(gfa::String, query_file::String)
    # Read the query ids from the file 
    query_ids = OrderedSet{String}()
    for line in eachline(query_file)
        push!(query_ids, strip(line))
    end
    queries = read_queries(gfa, query_ids)
    return queries, query_ids
end

function writeResults(ca::Vector{Int32}, color::Color, query_ids::OrderedSet{String}, out_file::String, size_map::Dict{Int32, Int32})
    h = open(out_file, "w+")
    prev_ori = Origin(-1,-1)
    aln_start = 1
    genome_loc = 1
    q_count = 1

    j = 0
    for i in 1:length(color.len) 
        i == 1 && println(color.origin[i])
        j +=1
        if ca[i] < 0 
            i > 1 && println(h, query_ids[q_count], "\t", color.origin[j-1].id, "\t", aln_start, "\t", genome_loc, "\t", color.len[j-1])
            q_count +=1
            genome_loc = 0
            prev_ori = Origin(-1,-1)

        elseif prev_ori.id != color.origin[i].id || prev_ori.pos != color.origin[i].pos 
            if genome_loc > 1 && color.len[i-1] > 0
                println(h, query_ids[q_count], "\t", color.origin[i-1].id, "\t", aln_start, "\t", genome_loc+(color.k_size - 1), "\t", color.len[i-1])
            end
            node_size = size_map[ca[i]]
            aln_start = copy(genome_loc)
            genome_loc +=  node_size - (color.k_size - 1)  # genome_loc == 1 ? node_size :
            prev_ori = color.origin[i]
        else 
            node_size = size_map[ca[i]]
            genome_loc +=  node_size - (color.k_size - 1)
        end 
    end
    # Do the last one
    println(h, query_ids[q_count-1], "\t", color.origin[j-1].id, "\t", aln_start, "\t", genome_loc, "\t", color.len[j-1])

end