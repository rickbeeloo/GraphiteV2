using OrderedCollections

include("./graphite.jl")
include("./graph_io.jl")
include("./suffixArray.jl")

@time run("/home/codegodz/Downloads/china_nl_graph.cf_seq", "/home/codegodz/Downloads/china_nl_graph.cf_seg", "/home/codegodz/SuffixTestRepo/dutch_queries.txt", Int32(31), "output.txt")