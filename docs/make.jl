using DelaySSAToolkit
using Documenter

DocMeta.setdocmeta!(DelaySSAToolkit, :DocTestSetup, :(using DelaySSAToolkit); recursive=true)

makedocs(;
    modules=[DelaySSAToolkit],
    authors="Xiaoming Fu",
    repo="https://github.com/palmtree2013/DelaySSAToolkit.jl/blob/{commit}{path}#{line}",
    sitename="DelaySSAToolkit.jl",
    format=Documenter.HTML(;
        mathengine=Documenter.Writers.HTMLWriter.MathJax2(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palmtree2013.github.io/DelaySSAToolkit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/tutorials.md", 
            "tutorials/examples.md"
            ],
        "Theory" => "theory.md",
        "Algorithm" => "Algorithms.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/palmtree2013/DelaySSAToolkit.jl",
    devbranch="main",
)
