import Pkg
projectPath = dirname(@__FILE__)
Pkg.activate(projectPath)
Pkg.instantiate()
