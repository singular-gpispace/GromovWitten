# # make file only for local make of the document.
# this result in errors in Travis
#
# ## To build the documentation locally
#
# ### Running inside Julia REPL:
#  cd to docs folder `cd("path to docs\\MatLang\\docs")` and run the following command:
# ```
# include("makeLocal.jl")
# ```
#
# ### Running inside OS Terminal:
# cd to docs folder using OS terminal and run the following command (julia path should be added to OS path):
# ```
# julia --color=yes makeLocal.jl
# ```
# Modify the LOAD_PATH to include the source directory
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Import Pkg module
using Pkg

# Activate the project in the parent directory
Pkg.activate(joinpath(@__DIR__, ".."))

# Install Documenter package
Pkg.add("Documenter")

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
        "Quasimodular" =>"quasimodular.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(;
#     repo="github.com/juliamatlab/MatLang.git",
# )
