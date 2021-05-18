using Pkg
pkg"add SumOfSquares#master"

using Documenter
using Literate
using CondensedMatterSOS
using Test

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

include(joinpath(EXAMPLES_DIR, "run_examples.jl"))

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
end

makedocs(
    sitename = "CondensedMatterSOS",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Introduction" => "index.md",
        "Examples" => map(EXAMPLES) do jl_file
            # Need `string` as Documenter fails if `name` is a `SubString{String}`.
            name = string(split(jl_file, ".")[1])
            return name => "generated/$name.md"
        end
    ],
)

deploydocs(
    repo   = "github.com/blegat/CondensedMatterSOS.jl.git",
    push_preview = true,
)
