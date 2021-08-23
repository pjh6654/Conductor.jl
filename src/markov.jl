abstract type AbstractMarkov <: AbstractKinetics end

get_output(x::AbstractMarkov) = x.states
openstates(x::AbstractMarkov) = x.open

struct Markov{S<:AbstractSystem} <: AbstractMarkov
    sys::S
    states::Vector
    ss::Union{Nothing,Vector}
    open::Vector
end

function Markov(rn::ReactionSystem, open::Vector{Num};
    sys_type::Type{S}=ODESystem,
    ss::Union{Bool,Vector{Num}}=false) where S<:AbstractSystem

    sys = convert(sys_type,rn)
    if ss isa Bool
        ss = ss ? _steady_state(rn) : nothing
    end
    if !isnothing(ss)
        display(get_states(rn))
        sys = @set sys.defaults = Dict(get_states(rn) .=> ss)
    end
    Markov(sys,get_states(sys),ss,open)
end

#=
Steady states of fully connected Markov Chains
Markov Kinetics can be written as a matrix equation
                dx/dt = Qx  (1)
where Q is the transition matrix (Jacobian of the system)

The steady state problem is as follows:
                dx/dt = Qx = 0  (2)
However Markov Kinetics for ion channels have a extra constraints
                xᵢ ≧ 0  (3)
                ∑xᵢ = 1  (4)
Therefore we can create an augmented matrix A using equation (4)
                A = let x = Q; x[end,:] = 1; x end
This produces a new system
                Ax = b
where b is:
                b = let x = zeros(size(A,1)); x[end] = 1; x end
Finding the unique steady state is as simple as finding
                x' = inv(A)b
                x' = rref([A b]) <- simpler output for symbolics
=#

# Row Reduced Echelon Form for symbolic calculations
function rref!(A::Matrix{T}) where {T<:Num}
    nr, nc = size(A)
    i = j = 1
    while i <= nr && j <= nc
        d = A[i,j]
        for k = j:nc
            A[i,k] /= d
        end
        for k = 1:nr
            if k != i
                d = A[k,j]
                for l = j:nc
                    A[k,l] -= d*A[i,l]
                end
            end
        end
        i += 1
        j += 1
    end
    A
end

rref(A::Matrix{T}) where {T<:Num} = rref!(copy(A))

function _steady_state(sys::ODESystem)
    # WARNING: This produces large symbolic equations
    Q = calculate_jacobian(sys)
    Q[end,:] .= 1
    b = zeros(Num,size(Q,1))
    b[end] = 1
    rref([Q b])[:,end]  # Row Reduced Echelon Form
end
_steady_state(sys::ReactionSystem) = _steady_state(convert(ODESystem,sys))

function _conductance(::Type{ODESystem},gbar_val::T, markov::AbstractMarkov;
                      passive::Bool=false,
                      null_init::Bool = false, name::Symbol) where {T <: Real}

    inputs = Set{Num}()
    states = Set{Num}()
    eqs = Equation[]

    markov_states = Set{Num}(x for x in get_output(markov))
    @variables g(t)
    push!(states, g)
    params = @parameters gbar
    defaultmap = Pair[gbar => gbar_val]

    syms = value.(union(get_variables.(equations(markov))...))
    for i in syms
        if i ∉ markov_states
            isparameter(i) ? push!(params, i) : push!(inputs, i)
        end
    end
    union!(states, markov_states, inputs)
    push!(eqs, g ~ gbar * sum(openstates(markov)))
    append!(eqs, equations(markov))
    if null_init
        foreach(markov_states) do x; defaultmap = _merge(defaultmap, Dict(x=>0.0)) end
    else
        default = hassteadystate(markov) ? defaults(markov) : Dict(markov_states .=> 1//length(markov_states)) # Sum of states must equal 1
        defaultmap = _merge(defaultmap, default)
    end

    system = ODESystem(eqs, t, states, params;defaults = defaultmap, name = name)
    return (collect(inputs), params, system)
end
