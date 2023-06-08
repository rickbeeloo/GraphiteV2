using IntervalTrees

struct genomeLocation
    ref_id::Int32 
    size::Int32 
end

const treeType = IntervalTrees.IntervalBTree{Int32,  IntervalValue{Int32, genomeLocation}, 64}
const treeIntersectType = IntervalTrees.IntervalIntersectionIterator{typeof(IntervalTrees.true_cmp), Int32, IntervalValue{Int32, String}, 64}


function _crosses_flank_simple(query::Tuple{Int32, Int32}, match::IntervalValue{Int32, genomeLocation})

    # Always with respect to the match

    first(query) <= first(match) && last(query) > first(match) && last(query) > last(match) && return :beyond

    first(query) < first(match) && last(query) >= first(match) && last(query) <= last(match) && return :right_flank 

    first(query) >= first(match) && first(query) <= last(match) && last(query) > last(match) && return :left_flank    
       
    first(query) >= first(match) && last(query) <= last(match) && return :engulfed

    first(query) == first(match) && last(query) == last(match) && return :identical

    println("WEIRD: ", query, " with match: ", match )

    error("Uknown intersection type")
end


function _update_tree(tree::treeType, query::Tuple{Int32, Int32}, query_size::Int32, match::IntervalValue{Int32, genomeLocation}, intersect_type::Symbol)

    (intersect_type == :engulfed || intersect_type == :identical) &&  return (0,0)
    
    if intersect_type == :beyond 
       # println("Engulfed")
        delete!(tree, (first(match), last(match)) )

    elseif intersect_type == :right_flank 
        if query_size > value(match).size
         #   println("RF: query longer")
            delete!(tree, (first(match), last(match)))
            tree[ (last(query) + Int32(1), last(match)) ] = value(match)
        else 
         #   println("RF, query shorter")
            # If we already shortened on the left before this might become a too short interval 
            # check this 
            new_right_side = first(match) - Int32(1)
            query = new_right_side > query[1] ? (query[1], new_right_side) : (0,0)
        end

    elseif intersect_type == :left_flank 
        if query_size > value(match).size
         #   println("LF: query longer")
            delete!(tree, (first(match), last(match)) )
            tree[first(match), first(query)-1] = value(match)
        else 
        #    println("LF Query shorter")
            # We have to shorten our query 
            query = (last(match) + Int32(1), query[2])
        end
    else 
        error("Shouldn't come here")    
    end

    return query
end

function _process_intersections(tree::treeType, intersections::Vector{IntervalValue{Int32, genomeLocation}}, interval::Tuple{Int32, Int32}, identifier::Int32, size::Int32)
    
    intersection_types = [_crosses_flank_simple(interval, intersection) for intersection in intersections]
   # println(intersection_types)
    
    for (i, intersection) in enumerate(intersections)
        #println("Intersection: ", intersection, " with interval: ", interval)
        intersect_type = intersection_types[i] #_crosses_flank_simple(interval, intersection)
       # println("------> ", intersection_types)
        interval = _update_tree(tree, interval, size, intersection, intersect_type)
    end
    #println("Interval at the end: ", interval)
    # Maybe we completely stripped down the interval
    if interval[1] != 0 && interval[2] != 0
        tree[interval] = genomeLocation(identifier, size)
    end 

end 

function color!(tree::treeType, interval::Tuple{Int32, Int32}, identifier::Int32, size::Int32)
    intersections = collect(intersect(tree, interval))
    if length(intersections) > 0
        _process_intersections(tree, intersections, interval, identifier, size)
    else
        tree[interval] = genomeLocation(identifier, size)
    end
    
end



# function main() 
#     xs = IntervalMap{Int32, genomeLocation}()

#     xs[( Int32(1), Int32(5) )]  = genomeLocation(    Int32(10), Int32(1))
#     xs[( Int32(6),  Int32(10))] = genomeLocation(    Int32(1),   Int32(1))
#     xs[  Int32(301),  Int32(390)]  = genomeLocation( Int32(1), Int32(1))
#     xs[( Int32(400),  Int32(410))] =genomeLocation(  Int32(11), Int32(1))

#     # println("Tree before")
#     # println(xs)
#     # println() 

#     color!(xs, ( Int32(301), Int32(401)), Int32(2), Int32(6))

#     # #color!(xs, ( Int32(6000), Int32(7000)), Int32(60))

#     println(xs)
# -
# end


# main()