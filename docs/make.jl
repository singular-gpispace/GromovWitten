#=
To build the documentation, do this from within the GromovWitten.jl directory:

julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl

To run the tests again, it suffices to repeat the second command.
=#

# Import required packages
using GromovWitten
using Documenter

# Set up documentation test setup
DocMeta.setdocmeta!(GromovWitten, :DocTestSetup, :(using GromovWitten; recursive=true))

makedocs(;
    modules=[GromovWitten],
    clean = true,
    checkdocs = :none,
    doctest = true,
    #strict = true,
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
        "Quasimodular" =>[
            "quasimodular.md",
            "quasimodular_psi.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/GromovWitten.git",
    devbranch="main",
)
