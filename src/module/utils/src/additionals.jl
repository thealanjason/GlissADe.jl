#=
# Typechecking for parsing.
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#

"""
            _check_valid_int(s)
Checks if the given variable is of type-constant `INT_TYPE`. 
## Arguments 
- s::AnyDataType - Any Variable 
"""
function _check_valid_int(s)
    try
        parse(INT_TYPE[], s)
        return true
    catch
        return false
    end
end

"""
            _check_valid_face(line)
Checks if the given string `line` represents a valid face. 
## Arguments
- line::String - line which needs to be checked 
"""
function _check_valid_face(line)
    parts = split(line, '(', limit = 2)
    length(parts) != 2 && return false
    @inbounds n_str = strip(parts[1])
    !_check_valid_int(n_str) && return false
    n = parse(INT_TYPE[], n_str)
    @inbounds face_str = strip(parts[2], [' ', ')'])
    face_elements = split(face_str)
    length(face_elements) != n && return false
    !all(_check_valid_int, face_elements) && return false
    return true
end
