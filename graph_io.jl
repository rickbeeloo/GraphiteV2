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

    # For later to speed up testing 
    q_paths = open("query_paths.txt", "r")

    queries = Vector{Vector{Int32}}()
    query_ids_file_order = OrderedSet{String}()

    for line in eachline(q_paths)
        identifier, path = split(line, "\t")
        push!(query_ids_file_order, identifier)
        push!(queries, path_numbers)
    end

    # found = 0
    # queries = Vector{Vector{Int32}}()
    # query_ids_file_order = OrderedSet{String}() # Might find in different order in the file than query_ids
    # for (i, line) in enumerate(eachline(f))
    #     identifier, path = split(line, "\t")
    #     if identifier in query_ids
    #         write(q_paths_out, line)
    #         found +=1
    #         path_numbers = parse_numbers(path)
    #         push!(query_ids_file_order, identifier)
    #         push!(queries, path_numbers)
    #         found == length(query_ids) && break
    #     end
    # end
    return queries, query_ids_file_order
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

function read_ids_from_file(f::String, first_n::Int64)
function read_ids_from_file(f::String, first_n::Int64)
    ids = OrderedSet{String}()
    for (i, line) in enumerate(eachline(f))
    for (i, line) in enumerate(eachline(f))
        push!(ids, strip(line))
        first_n > 0 && i == (first_n-1) && break 
        first_n > 0 && i == (first_n-1) && break 
    end
    return ids
end

function processGFA(gfa::String, query_file::String; first_n = -1)
function processGFA(gfa::String, query_file::String; first_n = -1)
    # Read the query ids from the file 
    query_ids = read_ids_from_file(query_file, first_n)
    # Get the first X if specified (handy for testing)
    query_ids = read_ids_from_file(query_file, first_n)
    # Get the first X if specified (handy for testing)
    queries, query_ids = read_queries(gfa, query_ids)
    return queries, query_ids
end

function writeResults(ca::Vector{Int32}, color::Color, query_ids::OrderedSet{String}, out_file::String, size_map::Dict{Int32, Int32})
    h = open(out_file, "w+")
    prev_ori =  color.origin[1]
    aln_start = 1
    genome_loc = 0
    q_count = 1
   
    for i in 1:length(color.len)-1
  
        if ca[i] > 0 
            node_size = size_map[ca[i]]
            genome_loc += node_size - color.k_size + 1
        end

        if ca[i+1] < 0 || prev_ori.id != color.origin[i+1].id || prev_ori.pos != color.origin[i+1].pos
            color.len[i] > 0 && println(h, query_ids[q_count], "\t", color.origin[i].id, "\t", aln_start, "\t", genome_loc+color.k_size-1, "\t", color.len[i])
            prev_ori = color.origin[i+1]
            aln_start = copy(genome_loc) + 1 
        end

        if ca[i+1] < 0
            q_count +=1 
            prev_ori = color.origin[i+1]
            genome_loc = 0
        end 

    end
        

end

