using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using GromovWitten
using Documenter
DocMeta.setdocmeta!(GromovWitten, :DocTestSetup, :(using GromovWitten); recursive=true)

makedocs(;
    modules=[GromovWitten],
    authors="Ali Traore",
    repo="https://github.com/singular-gpispace/GromovWitten.jl/blob/{commit}{path}#{line}",
    sitename="GromovWitten.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/GromovWitten.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/GromovWitten.git",
    devbranch="main",
)
