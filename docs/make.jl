using DelaySSAToolkit
using Documenter

DocMeta.setdocmeta!(
    DelaySSAToolkit, :DocTestSetup, :(using DelaySSAToolkit); recursive=true
)

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
            "tutorials/bursty.md",
            "tutorials/delay_degradation.md",
            "tutorials/heterogeneous_delay.md",
            "tutorials/delay_oscillator.md",
            "tutorials/stochastic_delay.md",
        ],
        "Algorithm" => [
            "algorithms/notations.md",
            "algorithms/delayrejection.md",
            "algorithms/delaydirect.md",
            "algorithms/delaymnrm.md",
        ],
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
    warnonly=true,
)

deploydocs(; repo="github.com/palmtree2013/DelaySSAToolkit.jl", devbranch="main")
