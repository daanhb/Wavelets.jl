# This file contains definitions of doubly infinite sequences

# generic interface
export shift

# interface for compact sequences
export zero_padding, CompactSequence
export support

# interface for periodic extensions
export PeriodicExtension
export period

# interface for symmetric extensions
export SymmetricExtension
export symmetric_extension_halfpoint_odd, symmetric_extension_halfpoint_even,
    symmetric_extension_wholepoint_odd, symmetric_extension_wholepoint_even

# interface for symmetric extensions
export CompactSymmetricExtension
export compact_symmetric_sequence_halfpoint_odd,
    compact_symmetric_sequence_halfpoint_even,
    compact_symmetric_sequence_wholepoint_odd,
    compact_symmetric_sequence_wholepoint_even


"""
A `DoublyInfiniteVector` is a sequence of infinite length that can be indexed
with any integer, including both positive and negative ones.

It is also referred to as a bi-infinite sequence or two-way infinite sequence.
"""
abstract type DoublyInfiniteVector{T}
end

const Sequence = DoublyInfiniteVector
const Filter = DoublyInfiniteVector

Base.eltype(::Type{DoublyInfiniteVector{T}}) where T = T
Base.eltype(::Type{V}) where V <: DoublyInfiniteVector = eltype(supertype(V))

Base.getindex(s::DoublyInfiniteVector, r::Range) = [s[i] for i in r]

# The parahermitian conjugate of a sequence is the conjugate of its reverse
parahermitian_conjugate(s::DoublyInfiniteVector) = conj(reverse(s))

adjoint(s::DoublyInfiniteVector) = conj(s)
