# FlexibilityAnalysis.jl

A JuMP extension for analyzing and quantifying the flexibility of complex systems.
Formerly known as [FlexJuMP.jl](https://github.com/pulsipher/FlexJuMP.jl).

## Status
`FlexibilityAnalysis.jl` is long out of date and is currently without a developer. It last worked with Julia `v1.1` and JuMP `v0.18` using Gurobi.jl `v.6` (wrapping Gurobi `0.7.x`). Please respond to [this issue](https://github.com/pulsipher/FlexibilityAnalysis.jl/issues/2) if you would like to help bring this up-to-date.

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/FlexibilityAnalysis.jl/stable) | [![Build Status](https://travis-ci.com/pulsipher/FlexibilityAnalysis.jl.svg?branch=v0.1.0)](https://travis-ci.com/pulsipher/FlexibilityAnalysis.jl) [![Build Status2](https://ci.appveyor.com/api/projects/status/github/pulsipher/FlexibilityAnalysis.jl?branch=master&svg=true)](https://ci.appveyor.com/project/pulsipher/FlexibilityAnalysis-jl) [![codecov.io](http://codecov.io/github/pulsipher/FlexibilityAnalysis.jl/coverage.svg?branch=master)](http://codecov.io/github/pulsipher/FlexibilityAnalysis.jl?branch=master) |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pulsipher.github.io/FlexibilityAnalysis.jl/dev) | [![Build Status](https://travis-ci.com/pulsipher/FlexibilityAnalysis.jl.svg?branch=master)](https://travis-ci.com/pulsipher/FlexibilityAnalysis.jl) [![Build Status2](https://ci.appveyor.com/api/projects/status/github/pulsipher/FlexibilityAnalysis.jl?branch=master&svg=true)](https://ci.appveyor.com/project/pulsipher/FlexibilityAnalysis-jl) [![codecov.io](http://codecov.io/github/pulsipher/FlexibilityAnalysis.jl/coverage.svg?branch=master)](http://codecov.io/github/pulsipher/FlexibilityAnalysis.jl?branch=master) |

Comments, suggestions and improvements are welcome and appreciated.

## License
`FlexibilityAnalysis` is licensed under the [MIT "Expat" license](./LICENSE).

## Installation
`FlexibilityAnalysis.jl` is a registered Julia package and can be installed in the usual manner.

```julia
using Pkg
Pkg.add("FlexibilityAnalysis")
```

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/FlexibilityAnalysis.jl/stable)

Please visit our [documentation pages](https://pulsipher.github.io/FlexibilityAnalysis.jl/stable) to learn all you need to know to get started and more. We provide
a quick start guide, an overview necessary background material, a detailed user guide, examples, and
an API library.

## Project Status
The package was tested (in 2019) against Julia `1.0`, `1.1`, and nightly on Linux, macOS, and Windows.

## Contributing
`FlexibilityAnalysis` needs a newer developer, see [this issue](https://github.com/pulsipher/FlexibilityAnalysis.jl/issues/2) to learn more.
