#=
The Cell struct is used to store local (face-only) geometrical data and 
field values (pressure, velocity, snow thickness)

Last Updated On: 11th January, 2025 20:53 UTC+5:30
=#

export Cell

"""
    mutable struct Cell{T,S,W} <: Abstract Cell 
Structure used for storing geometrical data. 

## DataTypes
- T - Type for numbers for geometrical data. Should be `Dual` when differentiating with geometry. 
- S - Type for integers. Defaults to `INT_TYPE`
- W - Type for numbers for state variables. Should be `Dual` when differentiation is performed. 

## Fields 
- idx::S - Face Index 
- center::Vector{T} - Coordinates of face centroid 
- vertices::Vector{Vector{T}} - Coordinates of the vertices forming the face. 
- edge_centers::Vector{Vector{T}} - Coordinates of edge centers of all edges of a face 
- edge_lengths::Vector{T} - Edge lengths of all edges of a face 
- normal::Vector{T} - Surface normal at the centroid of a face 
- area::T - Area of a face 
- edge_binormals::Vector{Vector{T}} - Binormals orthogonal to edges and pointing outwards for each edge of a face 
- transform::Vector{Matrix{T}} - Direction Cosine Matrix for transforming variable of current face for each edge. 
- transform2::Vector{Matrix{T}} - Direction Cosine Matrix for transforming variable of neighbouring face for each edge. 
- neighbours::Vector{S} - List of neighbouring faces. 
- h::W - Thickness at this cell 
- vel::Vector{W} - Velocity (in global coords) in this cell 
- pb::W - Basal Pressure at this cell 
"""
@with_kw mutable struct Cell{T, S, W} <: AbstractCell
    idx::S # index of the cell
    center::Vector{T} # coords of the centroid of the cell
    vertices::Vector{Vector{T}} # coords of the vertices of the cell -> A View of the global points array
    edge_centers::Vector{Vector{T}} # coords of all the edge centers
    edge_lengths::Vector{T} # Edge lengths of all edges of the cell
    normal::Vector{T} # Surface Normal
    area::T # Area of the cell
    edge_binormals::Vector{Vector{T}} # Edge Binormals [Central interpolation hardcoded]
    transform::Vector{Matrix{T}}
    transform2::Vector{Matrix{T}}
    neighbours::Vector{S} # Neighbours of the current cell
    h::W # Thickness
    vel::Vector{W} # Velocity
    pb::W # Basal Pressure
end
