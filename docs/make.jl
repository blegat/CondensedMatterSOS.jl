using Documenter
using Literate
using CondensedMatterSOS
using Test

const _EXAMPLE_DIR = joinpath(@__DIR__, "src", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

function link_example(content)
    edit_url = match(r"EditURL = \"(.+?)\"", content)[1]
    footer = match(r"^(---\n\n\*This page was generated using)"m, content)[1]
    content = replace(
        content, footer => "[View this file on Github]($(edit_url)).\n\n" * footer
    )
    return content
end

for file in readdir(_EXAMPLE_DIR)
    if !endswith(file, ".jl")
        continue
    end
    filename = joinpath(_EXAMPLE_DIR, file)
    # `include` the file to test it before `#src` lines are removed. It is in a
    # testset to isolate local variables between files.
    @testset "$(file)" begin
        include(filename)
    end
    Literate.markdown(filename, OUTPUT_DIR; documenter = true)
    Literate.markdown(filename, OUTPUT_DIR)
end

makedocs(
    sitename = "CondensedMatterSOS",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Introduction" => "index.md",
        "Examples" => map(
            file -> joinpath(OUTPUT_DIR, file),
            filter(
                file -> endswith(file, ".md"),
                sort(readdir(OUTPUT_DIR)),
            )
        )
    ],
)

deploydocs(
    repo   = "github.com/blegat/CondensedMatterSOS.jl.git",
)
