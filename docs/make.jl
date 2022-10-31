using tropicalfeynman
using Documenter

DocMeta.setdocmeta!(tropicalfeynman, :DocTestSetup, :(using tropicalfeynman); recursive=true)

makedocs(;
    modules=[tropicalfeynman],
    authors="Ali Traore",
    repo="https://github.com/blociss/tropicalfeynman.jl/blob/{commit}{path}#{line}",
    sitename="tropicalfeynman.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://blociss.github.io/tropicalfeynman.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/blociss/tropicalfeynman.jl",
    devbranch="main",
)
