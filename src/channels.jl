mutable struct AuxConversion
    params::Vector{Num}
    eqs::Vector{Equation}
end

# Conductance types (conductance as in "g")
abstract type AbstractConductance end

isbuilt(x::AbstractConductance) = x.sys !== nothing

struct IonChannel <: AbstractConductance
    gbar::SpecificConductance # scaling term - maximal conductance per unit area
    conducts::DataType # ion permeability
    inputs::Vector{Num} # cell states dependencies (input args to kinetics); we can infer this
    params::Vector{Num}
    kinetics::Union{AbstractKinetics,Vector{<:AbstractKinetics}} # kinetics equations
    sys::Union{AbstractSystem, Nothing} # symbolic system
end

# General purpose constructor
function IonChannel(conducts::Type{I},
                    gate_vars::Union{K,Vector{<:K}},
                    max_g::SpecificConductance = 0mS/cm^2;
                    name::Symbol,kwargs...) where {I <: Ion, K<:AbstractKinetics}
    # TODO: Generalize to other possible units (e.g. S/F)
    check_kinetics(gate_vars)
    sys_type = systemtype(gate_vars)

    gbar_val = ustrip(Float64, mS/cm^2, max_g)
    (inputs, params, system) = _conductance(sys_type,gbar_val, gate_vars; name = name, kwargs...)
    return IonChannel(max_g, conducts, inputs, params, gate_vars, system)
end

function (chan::IonChannel)(newgbar::SpecificConductance)
    newchan = deepcopy(chan)
    @set! newchan.gbar = newgbar
    gbar_val = ustrip(Float64, mS/cm^2, newgbar)
    gsym = length(newchan.kinetics) > 0 ?
           value(first(@parameters gbar)) :
           value(first(@parameters g))

    mapping = Dict([gsym => gbar_val])
    new_defaults = merge(defaults(newchan.sys), mapping)
    @set! newchan.sys.defaults = new_defaults
    return newchan
end

# Alias for ion channel with static conductance
function PassiveChannel(conducts::Type{I}, max_g::SpecificConductance = 0mS/cm^2;
                        sys_type::Type{<:AbstractSystem} = ODESystem,
                        name::Symbol = Base.gensym(:Leak)) where {I <: Ion}
    return IonChannel(conducts, Gate{sys_type}[], max_g; name = name, passive = true)
end

struct SynapticChannel <: AbstractConductance
    gbar::ElectricalConductance
    conducts::DataType
    reversal::Num
    inputs::Vector{Num}
    params::Vector{Num}
    kinetics::Vector{<:AbstractGate}
    sys::Union{ODESystem, Nothing}
end

function SynapticChannel(conducts::Type{I}, gate_vars::Vector{<:AbstractGate},
                         reversal::Num, max_g::ElectricalConductance = 0mS;
                         name::Symbol, kwargs...) where {I <: Ion}

    check_kinetics(gate_vars)
    sys_type = systemtype(gate_vars)

    gbar_val = ustrip(Float64, mS, max_g)
    (inputs, params, system) = _conductance(sys_type, gbar_val, gate_vars; null_init = true, name = name, kwargs...)
    return SynapticChannel(max_g, conducts, reversal, inputs, params, gate_vars, system)
end

function GapJunction(conducts::Type{I}, reversal::Num, max_g::ElectricalConductance = 0mS;
                     passive::Bool = false,
                     sys_type::Type{<:AbstractSystem} = ODESystem,
                     name::Symbol) where {I <: Ion}
    SynapticChannel(conducts, Gate{sys_type}[], reversal, max_g, name; passive = true)
end

function (chan::SynapticChannel)(newgbar::ElectricalConductance)
    newchan = deepcopy(chan)
    @set! newchan.gbar = newgbar
    gbar_val = ustrip(Float64, mS, newgbar)
    gsym = length(newchan.kinetics) > 0 ?
           value(first(@parameters gbar)) :
           value(first(@parameters g))

    mapping = Dict([gsym => gbar_val])
    new_defaults = merge(defaults(newchan.sys), mapping)
    @set! newchan.sys.defaults = new_defaults
    return newchan
end
