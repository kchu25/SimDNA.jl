using SimDNA
using Documenter

DocMeta.setdocmeta!(SimDNA, :DocTestSetup, :(using SimDNA); recursive=true)

makedocs(;
    modules=[SimDNA],
    authors="Shane Kuei Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/SimDNA.jl/blob/{commit}{path}#{line}",
    sitename="SimDNA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/SimDNA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/SimDNA.jl",
    devbranch="main",
)
