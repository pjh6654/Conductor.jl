# Return ODESystem pretty printing for our wrapper types
# FIXME: Large output from generated steady states for markov kinetics
#Base.show(io::IO, ::MIME"text/plain", x::IonChannel,::Type{<:AbstractGate}=eltype(x.kinetics)) = Base.display(isbuilt(x) ? x.sys : x)
