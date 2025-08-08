#=
Type-checking functions for internal use. 

Last Updated On: 11th January, 2025 20:34 UTC+5:30
=#


"""
            is_valid_integer(s)
Checks if the given variable is of type-constant `INT_TYPE`. 
## Arguments 
- s::AnyDataType - Any Variable 
"""
function is_valid_integer(s)
    global INT_TYPE
    try
        parse(INT_TYPE[], s)
        return true
    catch
        return false
    end
end

"""
            is_valid_face_line(line)
Checks if the given string `line` represents a valid face. 
## Arguments
- line::String - line which needs to be checked 
"""
function is_valid_face_line(line)
    parts = split(line, '(', limit = 2)
    length(parts) != 2 && return false
    @inbounds n_str = strip(parts[1])
    !is_valid_integer(n_str) && return false
    n = parse(INT_TYPE[], n_str)
    @inbounds face_str = strip(parts[2], [' ', ')'])
    face_elements = split(face_str)
    length(face_elements) != n && return false
    !all(is_valid_integer, face_elements) && return false
    return true
end
