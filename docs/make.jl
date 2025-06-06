using AstroCoords
using Documenter

DocMeta.setdocmeta!(AstroCoords, :DocTestSetup, :(using AstroCoords); recursive=true)

makedocs(;
    modules=[AstroCoords],
    authors="Alexander Belotserkovtsev <abelotserkovtsev@college.harvard.edu>, Anna Mikhailova ann.mihaylova0306@gmail.com",
    sitename="AstroCoords.jl",
    format=Documenter.HTML(;
        canonical="https://Alarion239.github.io/AstroCoords.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Alarion239/AstroCoords.jl",
    devbranch="main",
)
