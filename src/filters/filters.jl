module Filters

import Base: +, -, *, /

export zero_padding

include("sequences.jl")
include("extensions.jl")
include("operations.jl")
include("filterbank.jl")

end # module
