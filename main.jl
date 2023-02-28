using OrderedCollections
using ArgParse

include("./graphite.jl")
include("./graph_io.jl")
include("./suffixArray.jl")

# Function to parse command line options
function parse_commandline()
    printstyled("
   ___                 _     _ _       
  / _ \\_ __ __ _ _ __ | |__ (_) |_ ___ 
 / /_\\/ '__/ _` | '_ \\| '_ \\| | __/ _ \\
/ /_\\| | | (_| | |_) | | | | | ||  __/
\\____/|_|  \\__,_| .__/|_| |_|_|\\__\\___|
                |_|                    
                 
GRAPH MEM SEARCHER 2023
", color = :blue)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--query_file", "-q"
            help = "Path to file with query ids"
            required = true 
            arg_type = String
        "--gfa", "-g"
            help = "Path to gfa file"
            required = true
            arg_type = String
        "--seq", "-s"
            help = "Path to seq file"
            required = true
            arg_type = String
        "--k_size", "-k"
            help = "K-mer sized used to contract the graph"
            required = true
            arg_type = Int32
        "--out", "-o"
            help = "Output file path"
            required = true
            arg_type = String
    end
    return parse_args(s)
end

function main()  
    parsed_args = parse_commandline()
    println("Running aligner with the following arguments:")
    for (arg, val) in parsed_args
        println("  $arg  =>  $val")
    end
    run(parsed_args["gfa"], parsed_args["seq"], parsed_args["query_file"], Int32(parsed_args["k_size"]), parsed_args["out"])
end 

main()

#@time run("/home/codegodz/Downloads/china_nl_graph.cf_seq", "/home/codegodz/Downloads/china_nl_graph.cf_seg", "/home/codegodz/SuffixTestRepo/dutch_queries.txt", Int32(31), "output.txt")