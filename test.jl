function concat_with_seperator(vectors::Vector{Vector{Int32}},identifiers::Vector{Int32})
    length(vectors) == length(identifiers) || error("Vectors != identifiers")
    # Probably a cleaner way to do this :) Like using map 
    # to flip the node ids in slices and copying the slices to the output
    total_size = sum(map(length, vectors))
    concat_arr = zeros(Int32, total_size + length(vectors)) 
    vect_id = -1 * length(vectors) # to have it decending for sorting
    # Concat with seperator + store sign in most significant bit
    i = 1
    identifier_i = 1
    # Sometimes we might want to align back to the other queries just
    # not to itself. However the relationship between a sequence and it's suffixes 
    # is not directly evident as there are multiple suffixes along the SA that 
    # belong to the same query, likely all scattered. 
    query_offsets = Dict{Int32, UnitRange{Int32}}()
    last_vect_id = Int32(1)
    last_start = Int32(1)

    @inbounds for v in vectors
        for node_id in v 
            concat_arr[i] = node_id
            i += 1
        end 
        # Get the vector id offsets and store them in the query offsets dict 
        query_offsets[identifiers[identifier_i]] = last_start:i-1
        concat_arr[i] = vect_id 
        vect_id += 1
        i += 1
        last_vect_id = vect_id
        last_start = i
        identifier_i += 1
    end
    return concat_arr, query_offsets
end

concat_with_seperator( [Int32[1,2,3,4,5], Int32[100,200,100], Int32[10000,10000] ], [Int32(1), Int32(5), Int32(9)])