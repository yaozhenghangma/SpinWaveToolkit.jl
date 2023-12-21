using MagnonRenormalization
using Documenter

DocMeta.setdocmeta!(MagnonRenormalization, :DocTestSetup, :(using MagnonRenormalization); recursive=true)

makedocs(;
    modules=[MagnonRenormalization],
    authors="YaozhenghangMa <Yaozhenghang.Ma@gmail.com> and contributors",
    repo="https://github.com/YaozhenghangMa/MagnonRenormalization.jl/blob/{commit}{path}#{line}",
    sitename="MagnonRenormalization.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://YaozhenghangMa.github.io/MagnonRenormalization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YaozhenghangMa/MagnonRenormalization.jl",
    devbranch="main",
)
