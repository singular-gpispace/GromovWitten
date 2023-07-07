using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using TropicalFeynman
using Documenter
DocMeta.setdocmeta!(TropicalFeynman, :DocTestSetup, :(using TropicalFeynman); recursive=true)

makedocs(;
    modules=[TropicalFeynman],
    authors="Ali Traore",
    repo="https://github.com/singular-gpispace/TropicalFeynman.jl/blob/{commit}{path}#{line}",
    sitename="TropicalFeynman.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/TropicalFeynman.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/TropicalFeynman.git",
    devbranch="main",
)
