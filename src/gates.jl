abstract type AbstractGate <: AbstractKinetics end

hasexponent(x::AbstractGate) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false
get_output(x::AbstractGate) = x.state

struct Gate{S<:AbstractSystem} <: AbstractGate
    sys::S # symbolic system
    state::Num
    ss::Union{Nothing,Num}
    p::Float64 # optional exponent (defaults to 1)
end

abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::Type{SteadyStateTau}, name::Symbol, ss::Num, tau::Num, p::Float64;
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem) where S<:AbstractSystem
    sym = only(@variables $name(t))

    pars = Set{Num}()
    pars_defaults = Dict()
    for symbol in union(get_variables.([ss,tau])...)
        if isparameter(symbol)
            # Extra parameters are made unique to the gate ie. x => m₊x
            par = renamespace(name,symbol)
            # substitute in unique parameter
            ss, tau = substitute.([ss,tau], (symbol => par))
            # set parameter default
            haskey(defaults,symbol) ? push!(pars_defaults, par => defaults[symbol]) : nothing
            push!(pars,par)
        end
    end
    rn = ReactionSystem(
        [Reaction(ss/tau,nothing,[sym]),
        Reaction(1/tau,[sym],nothing)],
        t,[sym],pars;
        defaults=pars_defaults,
        name=name
    )
    sys = convert(sys_type,rn)
    Gate(sys,sym,ss,p)
end

function Gate(::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Float64;
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem) where S<:AbstractSystem
    sym = only(@variables $name(t))

    pars = Set{Num}()
    pars_defaults = Dict()
    for symbol in union(get_variables.([alpha,beta])...)
        if isparameter(symbol)
            # Extra parameters are made unique to the gate ie. x => m₊x
            par = renamespace(name,symbol)
            # substitute in unique parameter
            alpha, beta = substitute.([alpha,beta], (symbol => par))
            # set parameter default
            haskey(defaults,symbol) ? push!(pars_defaults, par => defaults[symbol]) : nothing
            push!(pars,par)
        end
    end
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    rn = ReactionSystem(
        [Reaction(alpha,nothing,[sym]),
        Reaction(alpha+beta,[sym],nothing)],
        t,[sym],pars;
        defaults=pars_defaults,
        name=name
    )
    sys = convert(sys_type,rn)
    Gate(sys,sym,ss,p)
end

Gate(t::Type{SteadyStateTau}, name::Symbol, alpha::Num, beta::Num, p::Real; kwargs...) =
Gate(t, name, alpha, beta, Float64(p); kwargs...)

Gate(t::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Real; kwargs...) =
Gate(t, name, alpha, beta, Float64(p); kwargs...)

# TODO: find a nicer way to do this
function Gate(::Type{SteadyStateTau}; p = one(Float64),
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem, kwargs...) where S<:AbstractSystem

    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau, :m, kwargs[:m∞], kwargs[:τₘ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau, :h, kwargs[:h∞], kwargs[:τₕ], p; defaults=defaults, sys_type=sys_type)
    else
        throw("invalid keywords")
    end
end

function Gate(::Type{AlphaBetaRates}; p = one(Float64),
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem, kwargs...) where S<:AbstractSystem

    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(AlphaBetaRates, :m, kwargs[:αₘ], kwargs[:βₘ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(AlphaBetaRates, :h, kwargs[:αₕ], kwargs[:βₕ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:αₙ, :βₙ], syms)
        Gate(AlphaBetaRates, :n, kwargs[:αₙ], kwargs[:βₙ], p; defaults=defaults, sys_type=sys_type)
    else
        throw("invalid keywords")
    end
end

function _conductance(::Type{ODESystem},gbar_val::T, gate_vars::Vector{<:AbstractGate};
                      passive::Bool = false, null_init::Bool = false,
                      name::Symbol) where {T <: Real}

    inputs = Set{Num}()
    states = Set{Num}()
    eqs = Equation[]
    # retrieve all variables present in the RHS of kinetics equations
    # g = total conductance (e.g. g(m,h) ~ ̄gm³h)
    if passive
        params = @parameters g
        defaultmap = Pair[g => gbar_val]
    else
        gates = Set{Num}(get_output(x) for x in gate_vars)
        @variables g(t)
        push!(states, g)
        params = @parameters gbar
        defaultmap = Pair[gbar => gbar_val]

        for i in gate_vars
            syms = value.(get_variables(only(equations(i))))
            for j in syms
                if j ∉ gates
                    isparameter(j) ? push!(params, j) : push!(inputs, j)
                end
            end
        end
        union!(states, gates, inputs)
        push!(eqs, g ~ gbar * prod(hasexponent(x) ? get_output(x)^x.p : get_output(x) for x in gate_vars))
        append!(eqs, only(equations(x)) for x in gate_vars)
        append!(defaultmap, get_output(x) => (hassteadystate(x) && !null_init) ? x.ss : 0.0 for x in gate_vars) # fallback to zero
        for x in gate_vars; defaultmap = _merge(defaultmap, defaults(x)) end
    end
    system = ODESystem(eqs, t, states, params; defaults = defaultmap, name = name)
    return (collect(inputs), params, system)
end
