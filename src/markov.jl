# Markov process
abstract type AbstractMarkovProcess end

output(x::AbstractMarkovProcess) = getfield(x, :output)
states(x::AbstractMarkovProcess) = getfield(x, :states)
transmatrix(x::AbstractMarkovProcess) = getfield(x, :transmatrix)
steadystate(x::AbstractMarkovProcess) = getfield(x, :steadystate)

struct MarkovProcess{T} <: AbstractMarkovProcess
    states::Vector{T}
    output::T
    transmatrix::Matrix{T}
    steadystate::Union{Nothing,Vector{T}}
end

function MarkovProcess(gates::GatingVariable...)
    states = markov_states(gates...)
    output = first(states)
    tmatrix = transmatrix(gates...)
    steadystate = markov_steadystate(gates...)
    states = map((x,y)->setdefault(x, y), states, steadystate)
    return MarkovProcess(states, output, tmatrix, steadystate)
end
MarkovProcess(gates::Vector{GatingVariable}) = MarkovProcess(gates...)

function markov_states(gates::GatingVariable...)
    name = ""
    indices = []
    for gate in gates
        name *= last(getmetadata(output(gate), Symbolics.VariableSource)) |> string
        push!(indices, Base.OneTo(gate.p+1))
    end
    name = Symbol(name)
    vars = eval(:(only(@variables $name[$(indices...)](t))))
    return ModelingToolkit.scalarize(vars) |> vec
end
markov_states(gates::Vector{GatingVariable}) = markov_states(gates...)

_kron_sum(A::AbstractMatrix, B::AbstractMatrix) = kron(A, oneunit(B)) + kron(oneunit(A), B)

function transmatrix(gate::GatingVariable)
    n = gate.p + 1
    A = zeros(Num, n, n)
    for i = 1:n
        A[i,i] = - (i-1)*gate.alpha - (n-i)*gate.beta
    end
    for i = 1:n-1
        A[i+1,i] = (n-i)*gate.beta
        A[i,i+1] = i*gate.alpha
    end
    return A
end
transmatrix(gates::GatingVariable...) = mapreduce(transmatrix, (x,y)->_kron_sum(x,y), reverse(gates))
transmatrix(gates::Vector{GatingVariable}) = transmatrix(gates...)

function markov_steadystate(gate::GatingVariable)
    n = gate.p
    alpha = gate.alpha
    beta = gate.beta

    steady = Vector{Num}(undef, n+1)
    denom = (alpha + beta)^n
    for i in 0:n
        num = binomial(n, i) * alpha^(n-i) * beta^i
        steady[i+1] = num/denom
    end
    return steady
end
markov_steadystate(gates::GatingVariable...) = mapreduce(markov_steadystate, (x,y)->kron(x,y), reverse(gates))
markov_steadystate(gates::Vector{GatingVariable}) = markov_steadystate(gates...)

# Ion Channel Constructors
function ConductanceSystem(g::Num, ion::IonSpecies, markov_process::MarkovProcess;
        gbar::Num, linearity::IVCurvature = Linear, extensions::Vector{ODESystem} = ODESystem[],
                           defaults = Dict(), name::Symbol = Base.gensym("Conductance"))

    eqs = Equation[]
    inputs = Set{Num}()
    markov_inputs = Set{Num}()
    embed_defaults = Dict()
    params = Set{Num}()

    isparameter(gbar) && push!(params, gbar)

    x = states(markov_process)
    Q = transmatrix(markov_process)
    foreach(state->push!(markov_inputs, state), x)
    _eqs = D.(x) .~ Q*x
    foreach(eq->get_variables!(inputs, eq), _eqs)
    append!(eqs, _eqs)

    for sym in inputs
        isparameter(sym) && push!(params, sym)
        hasdefault(sym) && push!(embed_defaults, sym => getdefault(sym))
    end
    # Remove parameters + generated states
    setdiff!(inputs, params, markov_inputs)
    push!(eqs, g ~ gbar * output(markov_process))
    sys = ODESystem(eqs, t, union(inputs, g, markov_inputs), params;
                    defaults = merge(embed_defaults, defaults), name = name)

    for ext in extensions
        sys = extend(sys, ext)
    end

    return ConductanceSystem(g, ion, GatingVariable[], inputs, sys, linearity, Q, name)
end

function IonChannel(ion::IonSpecies,
                    markov_process::MarkovProcess;
                    max_g::SpecificConductance = 0mS/cm^2,
                    extensions::Vector{ODESystem} = ODESystem[],
                    name::Symbol = Base.gensym("IonChannel"),
                    linearity::IVCurvature = Linear, defaults = Dict())
    if max_g isa SpecificConductance
        gbar_val = ustrip(Float64, mS/cm^2, max_g)
        @parameters gbar
        push!(defaults, gbar => gbar_val)
    else
        gbar = max_g
        if hasdefault(gbar)
            gbar_val = getdefault(gbar)
            if gbar_val isa SpecificConductance
                gbar_val = ustrip(Float64, mS/cm^2, gbar_val)
                gbar = setdefault(gbar, gbar_val)
            end
        end
    end
    @variables g(t)
    g = setmetadata(g, ConductorUnits, mS/cm^2) # TODO: rework with MTK native unit system
    ConductanceSystem(g, ion, markov_process;
                      gbar = gbar, name = name, defaults = defaults,
                      extensions = extensions, linearity = linearity)
end
