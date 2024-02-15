using SEMTor
using Documenter

DocMeta.setdocmeta!(SEMTor, :DocTestSetup, :(using SEMTor); recursive=true)

makedocs(;
    modules=[SEMTor],
    authors="Steffen Plunder <steffen.plunder@web.de> and contributors",
    sitename="SEMTor.jl",
    format=Documenter.HTML(;
        canonical="https://SteffenPL.github.io/SEMTor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SteffenPL/sEMTor.jl",
    devbranch="main",
)
