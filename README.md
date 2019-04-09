# FlexibilityAnalysis.jl

A JuMP extension for analyzing and quantifying the flexibility of complex systems.
Formerly known as [FlexJuMP.jl](https://github.com/pulsipher/FlexJuMP.jl).

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

The package is tested against Julia `1.0`, `1.1`, and nightly on Linux, macOS, and Windows. Julia `0.6`
is not supported by `FlexibilityAnalysis`, but the `julia-0.6` branch on [FlexJuMP.jl](https://github.com/pulsipher/FlexJuMP.jl)
can be used.

## Contributing
`FlexibilityAnalysis` is being actively developed and suggestions or other forms of contribution are encouraged.
There are many ways to contribute to this package:

- Suggest new/improved functionality
- Report an issue if you encounter some odd behavior, or if you have suggestions to improve the package.
- Contribute with code addressing some open issues, that add new functionality or that improve the performance.
- When contributing with code, add docstrings and comments, so others may understand the methods implemented.
- Contribute by updating and improving the documentation.
