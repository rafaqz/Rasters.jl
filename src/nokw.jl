# Often we need to distinguish between `nothing` and no user input,
# So that we can use automatic checks if the user didn't pass a keyword.
# NoKW is entirely for that. It should never be used from outside of the package.
struct NoKW end
const nokw = NoKW()
@inline isnokw(::NoKW) = true
@inline isnokw(_) = false
