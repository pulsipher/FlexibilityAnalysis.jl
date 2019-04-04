using Documenter, FlexibilityAnalysis, JuMP

makedocs(modules = [FlexibilityAnalysis],
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         sitename = "FlexibilityAnalysis.jl - A Framework for Flexibility Analysis",
         authors = "Joshua Pulsipher",
         pages = ["Home" => "index.md",
                  "Background" => "background.md",
                  "User Guide" => "guide.md",
                  "Examples" => "examples.md",
                  "Library" => "api.md"])

 deploydocs(
     repo   = "github.com/pulsipher/FlexibilityAnalysis.jl.git",
     target = "build",
     deps   = nothing,
     make   = nothing
 )
