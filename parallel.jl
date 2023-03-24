

function divide(queries::Vector{Vector{Int32}})
    threads = Threads.nthreads()
    total_size = sum(map(length, queries))
    println("Total size: ", total_size)
    println("Threads: ", threads)
    chunk_size_approx = floor(Int, total_size / threads)
    println("Approx: ", chunk_size_approx)
    chunk_size_approx > 0 || error("Bug in divide")
    println("Chunk size apprix: ", chunk_size_approx)

    # Get the indices to split on 
    chunk_starts = [1]

    # Get the query size and divide
    cur_len = 0
    for (i, query) in enumerate(queries)
        cur_len += length(query)
        if cur_len >= chunk_size_approx
            i + 1 > length(queries) && break
            push!(chunk_starts, i+1)
            cur_len = 0
        end 
    end
    push!(chunk_starts, length(queries)+1)

    # for ci in 1:length(chunk_starts)-1
    #     println(view(queries, chunk_starts[ci]:chunk_starts[ci+1]-1))
    # end
    
    return chunk_starts
end


# queries = [Int32[1,2,3,10,10,10,20,103,10], Int32[1,2], Int32[1,2,3,10,10,10,20], Int32[1,2]]
# println(divide(queries))