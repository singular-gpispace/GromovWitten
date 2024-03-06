# Modify the LOAD_PATH to include the source directory
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Import Pkg module
using Pkg

# Activate the project in the 'docs' folder
Pkg.activate(joinpath(@__DIR__, "..", "docs"))

# Update the project environment
Pkg.instantiate()

# Import required packages
using GromovWitten, Documenter

# Run doctests for GromovWitten
doctest(GromovWitten)
