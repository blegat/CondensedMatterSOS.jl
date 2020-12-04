const EXAMPLES = filter(ex -> endswith(ex, ".jl") && ex != "run_examples.jl" && ex != "all_examples.jl",
                        readdir(@__DIR__))
