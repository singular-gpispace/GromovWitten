using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "docs", "src"))  # Add docs/src/ to LOAD_PATH
using GromovWitten
using Documenter
DocMeta.setdocmeta!(GromovWitten, :DocTestSetup, :(using GromovWitten); recursive=true)

makedocs(;
    modules=[GromovWitten],
    clean = true,
    checkdocs = :none,
    doctest = true,
    strict = true,
    sitename="GromovWitten.jl",
    expandfirst = ["Overview.md"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/GromovWitten.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" =>[
            "index.md",
            "Installation.md", 
            "SmallExample.md",
        ],
        
        "Examples" =>[
            "Hurwitz.md",
            "Gromov.md",
            #"Feynman Integral/Feynman.md",
        ],
        "Functions" => "Overview.md",
        #"Quasimodular" =>"Feynman Integral/Feynman.md"
        "Quasimodular" =>"quasimodular.md"
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/GromovWitten.git",
    branch = "gh-pages",
    devbranch="dev",
    devurl = "dev",

)
