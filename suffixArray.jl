

const LIBSAIS = "libsais-2.7.1/libsais.so.2"


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

function create_suffix_array(in_vector::Vector{Int32}, free_space::Int32)
    out_vector = zeros(Int32, length(in_vector) + free_space)
    d = Dict(Iterators.map(reverse,pairs(sort(unique(in_vector)))))
    in_vector = Int32.(get.(Ref(d), in_vector, 0))
    n = length(in_vector)
    k = maximum(in_vector) +1 
    @ccall LIBSAIS.libsais_int(in_vector::Ptr{Int32}, out_vector::Ptr{Int32}, n::Int32, k::Int32, free_space::Int32)::Int32
    out_vector .+= 1
    return out_vector
end

function create_k_suffix_array(vectors::Vector{Vector{Int32}}, query_identifiers::Vector{Int32}, free_space::Int32)
    concat_array, query_offsets = concat_with_seperator(vectors, query_identifiers)
    suffix_array = create_suffix_array(concat_array, free_space)
    return concat_array, query_offsets, suffix_array
end

function build_lcp(sa::Vector{Int32}, V::Vector{Int32}, base::Integer=1)
    T = eltype(sa)
    pos = sa .+ T(1-base)
    n = length(pos)
    lcparr = similar(pos)
    rank = invperm(pos)
    h = 0
    for i in 1:n
        if rank[i] == 1
            continue
        end
        j = pos[rank[i]-1]
        maxh = n - max(i, j)
        while h <= maxh && V[i+h] == V[j+h]
            h += 1
        end
        lcparr[rank[i]] = h
        h = max(h-1, 0)
    end
    lcparr[1] = 0
    return lcparr
end

function locate_insert_point(sa::Vector{Int32}, concat_arr::Vector{Int32}, ref::AbstractVector{Int32})
    low = 1 
    high = length(sa)
    suffix_arr_len = length(sa)
    while low <= high 
        mid = low + ((high - low) >>> 0x01) 
        suffix_start = sa[mid]
        suffix = view(concat_arr, suffix_start:suffix_arr_len)
        if ref < suffix 
            high = mid - 1
        elseif ref > suffix
            low = mid + 1
        else 
            return Int32(mid)  # then they should be equal, exact match
        end
    end
    return Int32(low)
end

function inverse_perm_sa(sa::Vector{Int32})
    inv_sa_perm = similar(sa)
    for i in 1:length(sa)
        inv_sa_perm[sa[i]] = i
    end
    return inv_sa_perm
end