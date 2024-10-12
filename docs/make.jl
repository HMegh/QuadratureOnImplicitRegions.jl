using Documenter
using QuadratureOnImplicitRegions

makedocs(
    sitename = "QuadratureOnImplicitRegions",
    format = Documenter.HTML(),
    modules = [QuadratureOnImplicitRegions],
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Reference" => "reference.md"
    ]

)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/HMegh/QuadratureOnImplicitRegions.jl.git",
)