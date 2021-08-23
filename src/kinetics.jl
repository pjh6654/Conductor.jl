abstract type AbstractKinetics end

system(x::AbstractKinetics) = x.sys
states(x::AbstractKinetics) = union(get_output(x),get_states(system(x)))
equations(x::AbstractKinetics) = get_eqs(system(x))
defaults(x::AbstractKinetics) = defaults(system(x))
hassteadystate(x::AbstractKinetics) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false

systemtype(x::AbstractKinetics) = typeof(system(x))
systemtype(x::Type{K}) where K<:AbstractKinetics = fieldtype(K,:sys)
systemtype(x::Vector{K}) where K<:AbstractKinetics = systemtype(eltype(x))

check_kinetics(x::AbstractKinetics) = @assert systemtype(x) <: AbstractSystem "Kinetics must contain an AbstractSystem!"
function check_kinetics(x::Vector{<:AbstractKinetics})
    check_kinetics.(x)
    @assert ~isequal(systemtype(x), AbstractSystem) "Kinetics must contain the same AbstractSystem subtype!"
end
