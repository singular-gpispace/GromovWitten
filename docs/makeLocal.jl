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
using Pkg
pkg"activate .."
push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

#

using Documenter, GromovWitten

makedocs(;
    modules=[GromovWitten],
	format = Documenter.HTML(
        prettyurls = prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singular-gpispace/GromovWitten.jl/{commit}{path}#{line}",
    sitename="GromovWitten",
    authors="Ali Traore",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(;
#     repo="github.com/juliamatlab/MatLang.git",
# )
