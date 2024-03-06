
# Modify the LOAD_PATH to include the source directory
#push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
# Activate the project in the parent directory
#Pkg.activate(joinpath(@__DIR__, "..", "docs")) # Activate the 'docs' project

# Import Pkg module
using Pkg

# Activate the project in the parent directory
#Pkg.activate(joinpath(@__DIR__, ".."))

# Install Documenter package
#Pkg.add("Documenter")

# Import required packages
using GromovWitten
#using Documenter

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
        "Quasimodular" =>"quasimodular.md"
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/GromovWitten.git",
    devbranch="main",
)
