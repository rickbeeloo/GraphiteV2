using OrderedCollections
using Profile, FileIO
using LoopVectorization
using ProgressMeter

#include("./graph_io.jl")
include("./suffixArray.jl")

const MASK = Int32(1<<30) 

flipnode(n::Int32) = n ⊻ MASK
isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) 
get_original(n::Int32) = isflipped(n) ? flipnode(n) : flipnode(n)

struct Origin 
    id::Int32 
    pos::Int32 
end

struct Color 
    len::Vector{UInt32}
    origin::Vector{Origin}
    size_map::Dict{Int32, UInt32}
    k_size::Int32
end


function get_match_size(color::Color, ca::Vector{Int32}, match_start::Int32, match_size::Int32)
    # This would alloc cause of intermediate array
    #match_size_nt = sum(get.(Ref(color.size_map), view(ca, match_start:match_end), 0))
    match_size_nt = 0
    match_end = match_start+match_size-1
    for i in match_start:match_end
        match_size_nt += get(color.size_map, ca[i], 0)
    end
    # println("MS: ", match_size)
    match_size_nt = match_size_nt - ((match_size-1) * (color.k_size-1)) 
    return match_size_nt
end


function update_color!(color::Color, ref_id::Int32, match_start::Int32, match_size::Int32, ca::Vector{Int32})
    
    match_end = match_start+match_size-Int32(1)
        
    # Get info from pervious matches
    at_start = color.origin[match_start]
    at_end   = color.origin[match_end]

    if at_start.id > 0 && at_start.id == at_end.id && at_start.pos == at_end.pos
        return
    # Don't have to bother about single node matches as they can't be longer anyway
    elseif match_start ==  match_end && color.len[match_start] > 0
        return
    else
        match_size_nt =  get_match_size(color, ca, match_start, match_size) #sum(get.(Ref(color.size_map), view(ca, match_start:match_end), 0))
        #println("COLOR_UPDATE: adding size: ", match_size_nt)
        @inbounds @simd  for i in match_start:match_end
            if color.len[i] < match_size_nt 
                color.len[i]  = match_size_nt
                color.origin[i] = Origin(ref_id, match_start)
            end
        end
    end
end

function reverse_complement_ref!(ref::Vector{Int32})
    reverse!(ref)
    @inbounds for i in eachindex(ref)
        ref[i] = flipnode(ref[i])
    end
end

function convert_nodes!(in_vector::Vector{Int32})
    for i in eachindex(in_vector)
        in_vector[i] = convert_node(in_vector[i])
    end
end

# Check left and right from insert point to see if we have a match
function decide_on_seed(insert_point::Int32, ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32)
    # Check left for a match
    left_of_ip = insert_point > 1 ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point-Int32(1),  Int32(0))[1] : 0 
    left_of_ip > 0 && return insert_point-Int32(1), Int32(left_of_ip)

    # Check right for a match, no need to check if it's outside the bounds of the SA <= length(sa)
    right_of_ip = insert_point <= length(sa) ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point, Int32(0))[1] : 0
    right_of_ip > 0 && return insert_point, Int32(right_of_ip)

    # Neither actually matches our Ref, return 0,0 to move on to the next one
    return Int32(0), Int32(0)
end

function matches_till(ref::AbstractVector{Int32}, ref_start::Int32, ca::Vector{Int32}, q_start::Int32)
    (ref_start > length(ref) || q_start > length(ca)) && return 0
    smallest_n = min(length(ref)-ref_start+1, length(ca)-q_start+1)
    for i in 1:smallest_n
        if ref[ref_start + i - 1] != ca[q_start+i-1]
            return Int32(i - 1)
        end 
    end 
    return Int32(smallest_n)
end

function check_this_point(ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32, point::Int32, skip::Int32)
    # Given a point in the suffix array, compare the suffix to the Ref 
    ca_suffix_index = sa[point]
    ca_start = ca_suffix_index + skip
    #println("Checking suffix: ", view(ref, ref_start:length(ref)))
    #println("Checking suffix: ", view(ref, ref_start:length(ref)))
    ref_start = ref_start + skip
    match_size = matches_till(ref, ref_start, ca, ca_start) + skip
    return match_size, match_size - skip
end

function extend_from_point!(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, lcp::Vector{Int32}, point::Int32, forward::Bool, ref_start::Int32, match_size::Int32, ref_id::Int32, color::Color)
    
    moves = 0 
    checks = 0
    
    move_dir = forward ? 1 : -1
    lcp_dir  = forward ? 0 :  1
    i = point += move_dir

    println("Extending for: ", ref_id, " from position: ", ref_start)
    while i > 1 && i <= length(sa) && lcp[i + lcp_dir] > 0
        moves += 1
        # We can skip the LCP part when extending, note though we also have to 
        # check the previous match size so min(lcp valu, prev match size)
        start_check_from = Int32(min(lcp[i + lcp_dir], match_size))
        #println("Checking in SA at: ", i," with skip of: ", start_check_from )
        # Check the size of this match starting from +1 of the LCP value)
        #println("-- Suff: ", view(ref, ref_start:ref_start+20))
        match_size, actual_checks = check_this_point(ca, sa, ref, ref_start, Int32(i), start_check_from)
        println("Checking in SA at: ", i," with skip of: ", start_check_from, " matched: ", match_size, " checks done: ", actual_checks)
        checks += actual_checks
        update_color!(color, ref_id, sa[i], Int32(match_size), ca)
        i += move_dir     
        println("Condition checks")
        println(i <= length(sa))
        println(lcp[i + lcp_dir] > 0)   
    end
    println()
    return match_size, moves, checks
end

function align(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
   # println(ref_id, " at: ", ref[307])
    max_match_index = Int32(0)
    max_match_size = Int32(0)
    ref_start = Int32(1)

    # Keep track of some stats 
    bin_searches = 0
    nodes_updated = 0
    tot_moves = 0
    tot_checks = 0

    while ref_start <= length(ref)

        bin_searches += 1 
        println("Doing binary search again")
        
        insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
        max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)

        # If we have a match keep using the linked list to extend 
        if max_match_size > 0 
            while ref_start <= length(ref)
                println("Working at ref position: ", ref_start)
                
                # Check the match size at this point 
                max_match_size, _ = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size)-Int32(1)) # skip k-1
                                
                # If we don't have any match we don't have to check the flanks
                max_match_size == 0 && break 
                
                println("> Forward")
                updates, moves, checks = extend_from_point!(ca, sa, ref, lcp, max_match_index, true, ref_start, Int32(max_match_size), ref_id, color)
                nodes_updated += updates
                tot_moves += moves
                tot_checks += checks
                
                println("> Backward")
                updates, moves, checks =  extend_from_point!(ca, sa, ref, lcp, max_match_index, false, ref_start, Int32(max_match_size), ref_id, color)
                nodes_updated += updates
                tot_moves += moves
                tot_checks += checks

                # Move to next location in suffix array for the second around
                max_match_index = inv_perm_sa[sa[max_match_index]+1]
                ref_start += Int32(1)
            end 
        else 
            # No match at current point, move +1 to do a binary search again
            ref_start += Int32(1)
        end
    end

    println("Stats for: ", ref_id)
    println("Bin searches: ", bin_searches)
    println("Nodes updated: ", nodes_updated)
    println("Moves: ", tot_moves)
    println("Checks: ", tot_checks)
    println()

end

function align_forward_and_reverse(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    # First do the forward align 
    println("REF: ", ref_id, " with size: ", length(ref))
    println("> Forward")
    @time align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)

    # Flip the nodes and reverse to do the reverse alignment 
    reverse_complement_ref!(ref)
    println("> Reverse")
    @time align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)
    println()  
end

function run(gfa::String, seq_file::String, query_file::String, k_size::Int32, out_file::String; blacklist::String = "") 
    
   #blacklist_ids = !isempty(blacklist) ? read_ids_from_file(blacklist) : OrderedSet{String}()
    blacklist_ids = OrderedSet{String}()

    println("Reading queries")
    queries, query_ids = processGFA(gfa, query_file)
    
    println("Creating suffix array")
    ca, sa = create_k_suffix_array(queries, Int32(0))
    println("CA size: ", length(ca))
    println("Creating inv perm")
    inv_sa_perm = inverse_perm_sa(sa)
    println("Building LCP")
    lcp = build_lcp(sa, ca)

    println("Get node sizes")

    size_map = read_node_sizes(seq_file)
    #size_map = Dict( unique([(queries...)...]) .=> 50)
    #println(size_map)
    #size_map = Dict( unique([(queries...)...]) .=> 50)
    #println(size_map)
    len = zeros(Int32, length(ca))
    ori =  [Origin(-1,-1) for i in 1:length(ca)] # Vector{Origin}(undef, length(ca)) 
    color = Color(len, ori, size_map, k_size)

    limit = 1_800_000
    p = Progress(limit)
    println("Start aligning...")
    for (ref_id, line) in enumerate(eachline(gfa))
        if ref_id == 370
            println("Got the ref")
            identifier, path = split(line, "\t")
            ref_id == limit && break
            if !(identifier in query_ids) && !(identifier in blacklist_ids)
                path_numbers = parse_numbers(path)
                align_forward_and_reverse(Int32(ref_id), color, ca, sa, path_numbers, inv_sa_perm, lcp)
                next!(p)
            end
        end
    end

    #writeResults(ca, color, query_ids, out_file, size_map)
    save("test.jlprof",  Profile.retrieve()...)
end

#run("test", "test", "test", Int32(31), "test")

