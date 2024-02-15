using sEMTor
using Documenter

DocMeta.setdocmeta!(sEMTor, :DocTestSetup, :(using sEMTor); recursive=true)

makedocs(;
    modules=[sEMTor],
    authors="Steffen Plunder <steffen.plunder@web.de> and contributors",
    sitename="sEMTor.jl",
    format=Documenter.HTML(;
        canonical="https://SteffenPL.github.io/sEMTor.jl",
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
