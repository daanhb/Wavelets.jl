# This file contains definitions of sequences that arise from the extension
# of a finite vector.

# Utility function: return the first valid index of a generic vector
firstindex(v::AbstractVector) = first(eachindex(v))


"The supertype of doubly infinite sequences obtained by extending finite data."
abstract type ExtensionSequence{T} <: DoublyInfiniteVector{T}
end

"The length of the data vector that is stored in the sequence."
data_length(s::ExtensionSequence) = length(unsafe_data(s))

data(s::ExtensionSequence) = copy(unsafe_data(s))

"The leftmost index of the embedded data in the sequence."
leftindex(s::ExtensionSequence) = offset(s)
"The rightmost index of the embedded data in the sequence."
rightindex(s::ExtensionSequence) = offset(s) + data_length(s) - 1
"The support of the embedded data in the sequence."
data_support(s::ExtensionSequence) = leftindex(s):rightindex(s)

"The convex hull of the supports of the embedded data in two sequences."
function joint_support(s1::ExtensionSequence, s2::ExtensionSequence)
    l = min(leftindex(s1), leftindex(s2))
    r = max(rightindex(s1), rightindex(s2))
    l:r
end

"Shift the sequence `k` positions to the right."
shift(s::ExtensionSequence, k::Int) = similar(s, data(s), offset(s)+k)

# The extension sequence acts as a mutating view, elements of the underlying
# data can be altered.
Base.setindex!(s::ExtensionSequence, val, k::Int) =
    setindex!(unsafe_data(s), val, mapindex(s, k))

Base.widen(s::ExtensionSequence) = similar(s, widen(unsafe_data(s)), offset(s))

# Take element-wise conjugates
Base.conj(s::ExtensionSequence) = similar(s, conj(unsafe_data(s)), offset(s))

# We can do arithmetics generically by performing operations on the underlying
# data. Concrete sequences should implement `align` which aligns the data so
# that this is possible.
for op in (:+, :-)
    @eval function $op(s1::S, s2::S) where S <: ExtensionSequence
        a, b = align(s1, s2)
        similar(s1, $op(unsafe_data(a), unsafe_data(b)), offset(a))
    end
end

for op in (:+, :-, :*, :/)
    @eval ($op)(s::ExtensionSequence, x::Number) = similar(s, ($op)(unsafe_data(s),x), offset(s))
end

for op in (:+, :-, :*)
    @eval ($op)(x::Number, s::ExtensionSequence) = ($op)(s, x)
end

-(s::ExtensionSequence) = similar(s, -unsafe_data(s), offset(s))



################
## Sequences with compact support, or extension by zero padding
################

# For documentation see constructor below
struct CompactSequence{T} <: ExtensionSequence{T}
    # We store the non-zero coefficients as a regular vector
    v       ::  Vector{T}
    # offset stores the index at which the data starts, default is 1
    offset  ::  Int
end

const FIRFilter = CompactSequence

"""
    CompactSequence(v[, offset])

A `CompactSequence` is a compactly supported sequence, i.e., a sequence with
a finite number of nonzero elements. Its support is an interval ``[i,j]`` where
both `i` and `j` are integers and `j >= i`. The sequence elements are zero
outside that interval.

The `CompactSequence` stores a vector of length `j-i+1`. It can be thought of
as an extension of this vector to doubly infinite sequences by zero padding.
The optional `offset` in the constructor is the left endpoint `i` of the support
of the sequence. By default, it equals the first valid index of `v`.
"""
CompactSequence(v::Vector) = CompactSequence(v, 1)
CompactSequence(v::AbstractVector, offset = firstindex(v)) =
    CompactSequence(collect(v), offset)

"Extend the given vector to a doubly infinite sequency by zero padding."
zero_padding(v::AbstractVector, optional...) = CompactSequence(v, optional...)

unsafe_data(s::CompactSequence) = s.v
offset(s::CompactSequence) = s.offset

Base.similar(s::CompactSequence, v::AbstractVector, offset::Int) =
    CompactSequence(v, offset)

"A range of indices that includes all non-zero elements of the sequence."
support(s::CompactSequence) = data_support(s)

# For internal use:
# - map the index k of a sequence into an index l of the data
mapindex(s::CompactSequence, k) = k - offset(s) + 1
# - and vice-versa
imapindex(s::CompactSequence, l) = l + offset(s) - 1

# We override getindex to return zero outside our embedded vector.
Base.getindex(s::CompactSequence, k::Int) =
    k < leftindex(s) || k > rightindex(s) ? zero(eltype(s)) : getindex(unsafe_data(s), mapindex(s, k))

# Reverse the sequence in time
Base.reverse(s::CompactSequence) = CompactSequence(flipdim(unsafe_data(s), 1), -rightindex(s))

Base.ctranspose(s::CompactSequence) = CompactSequence(conj(flipdim(unsafe_data(s),1)), -rightindex(s))

# Return two sequences that have the same support, by zero padding
function align(s1::CompactSequence, s2::CompactSequence)
    if support(s1) == support(s2)
        # they already align, do nothing
        s1, s2
    else
        #they do not align, make new sequences
        supp = joint_support(s1, s2)
        CompactSequence(s1[supp], first(supp)), CompactSequence(s2[supp], first(supp))
    end
end


################
## Periodic extension
################

# For documentation see constructor below
struct PeriodicExtension{T} <: ExtensionSequence{T}
    v       ::  Vector{T}
    offset  ::  Int
end

"""
    PeriodicExtension(v[, offset])

Extend a given vector `v` using periodicity to a doubly infinite sequence.

The `offset` is optional and by default equals the first index of `v`.
"""
PeriodicExtension(v::Vector) = PeriodicExtension(v, 1)
PeriodicExtension(v::AbstractVector, offset = firstindex(v)) =
    PeriodicExtension(collect(v), offset)

periodic_extension(v::AbstractVector, optional...) = PeriodicExtension(v, optional...)

unsafe_data(s::PeriodicExtension) = s.v
offset(s::PeriodicExtension) = s.offset

period(s::PeriodicExtension) = data_length(s)

Base.similar(s::PeriodicExtension, v::AbstractVector, offset::Int) =
    PeriodicExtension(v, offset)

# For internal use:
# - map the index k of a sequence into an index l of the data
mapindex(s::PeriodicExtension, k) = mod(k - offset(s), period(s)) + 1
# - and vice-versa
imapindex(s::PeriodicExtension, l) = l + offset(s) - 1

Base.getindex(s::PeriodicExtension, k::Int) = getindex(unsafe_data(s), mapindex(s, k))

Base.reverse(s::PeriodicExtension) =
    PeriodicExtension(flipdim(unsafe_data(s), 1), -rightindex(s))

# Return two periodic sequences that are aligned (have the same length and offset)
function align(s1::PeriodicExtension, s2::PeriodicExtension)
    @assert period(s1) == period(s2)
    if offset(s1) == offset(s2)
        s1, s2
    else
        s1, PeriodicExtension(circshift(unsafe_data(s2), offset(s2)-offset(s1)), offset(s1))
    end
end


################
## Symmetric extensions
################


"""
A `SymmetricExtension` extends a vector `v` symmetrically to a doubly infinite
sequence.

The symmetry around each of the endpoints can be whole-point (the endpoint is
not repeated) or half-point (the endpoint is repeated). The symmetry can also be
even (symmetric) or odd (anti-symmetric).

Parameters:
- PL : either :wp (whole point) or :hp (half point) near left endpoint
- PR : eith :wp or :hp for right endpoint
- SL : either :odd or :even symmetry at left endpoint
- SR : either :odd or :even at right endpoint
"""
struct SymmetricExtension{T,PL,PR,SL,SR} <: ExtensionSequence{T}
    v       ::  Vector{T}
    offset  ::  Int

    SymmetricExtension{T,PL,PR,SL,SR}(v::AbstractVector{T}, offset = first(eachindex(v))) where {T,PL,PR,SL,SR} =
        new(collect(v), offset)
    SymmetricExtension{T,PL,PR,SL,SR}(v::Vector{T}) where {T,PL,PR,SL,SR} =
        new(v, first(eachindex(v)))
    SymmetricExtension{T,PL,PR,SL,SR}(v::Vector{T}, offset::Int) where {T,PL,PR,SL,SR} =
        new(v, offset)
end

unsafe_data(s::SymmetricExtension) = s.v
offset(s::SymmetricExtension) = s.offset

Base.similar(s::SymmetricExtension{T,PL,PR,SL,SR}, v::AbstractVector{S}, offset::Int) where {S,T,PL,PR,SL,SR} =
    SymmetricExtension{S,PL,PR,SL,SR}(v, offset)

# Provide four of the sixteen combinations for convenience. The other combinations
# can be constructed by explicitly calling the full constructor.
symmetric_extension_wholepoint_even(v::AbstractVector{T}) where T =
    SymmetricExtension{T,:wp,:wp,:even,:even}(v)

symmetric_extension_halfpoint_even(v::AbstractVector{T}) where T =
    SymmetricExtension{T,:hp,:hp,:even,:even}(v)

symmetric_extension_wholepoint_odd(v::AbstractVector{T}) where T =
    SymmetricExtension{T,:wp,:wp,:odd,:odd}(v)

symmetric_extension_halfpoint_odd(v::AbstractVector{T}) where T =
    SymmetricExtension{T,:hp,:hp,:odd,:odd}(v)

# Flip index k around center c using whole point symmetry:
#    c + (c - k) = 2*c-k
symmetric_flip_right(c, k, ::Val{:wp}) = 2*c-k
symmetric_flip_left(c, k, ::Val{:wp}) = 2*c-k

# Flip index k around center c using half point symmetry:
#    c+1/2 + (c+1/2 - k) = 2*c-k+1
symmetric_flip_right(c, k, ::Val{:hp}) = 2*c-k+1
#    c-1/2 + (c-1/2 - k) = 2*c-k-1
symmetric_flip_left(c, k, ::Val{:hp}) = 2*c-k-1

# Compute the index by mapping any index outside the range of the embedded vector
# to an index that is closer to the interval (using symmetry) and repeat.
# The recursion ends when the index lands inside the interval, which is hopefully
# quickly.
function mapindex(s::SymmetricExtension{T,PL,PR,SL,SR}, k) where {T,PL,PR,SL,SR}
    if k > rightindex(s)
        # We are to the right of the interval: use symmetry wrt right endpoint
        mapindex(s, symmetric_flip_right(rightindex(s), k, Val{PR}()))
    elseif k < leftindex(s)
        # We are to the left of the interval
        mapindex(s, symmetric_flip_left(leftindex(s), k, Val{PL}()))
    else
        k - offset(s) + 1
    end
end

parity_sign(::Val{:odd}) = -1
parity_sign(::Val{:even}) = 1

# For getindex we have to use the same logic as mapindex, but now we also have
# to take the parity (odd/even) of the symmetries into account in order to
# multiply by +/- 1.
function Base.getindex(s::SymmetricExtension{T,PL,PR,SL,SR}, k) where {T,PL,PR,SL,SR}
    if k > rightindex(s)
        # We are to the right of the interval: use symmetry wrt right endpoint
        parity_sign(Val{SR}()) * getindex(s, symmetric_flip_right(rightindex(s), k, Val{PR}()))
    elseif k < leftindex(s)
        # We are to the left of the interval
        parity_sign(Val{SL}()) * getindex(s, symmetric_flip_left(leftindex(s), k, Val{PL}()))
    else
        getindex(unsafe_data(s), k-offset(s)+1)
    end
end

function Base.reverse(s::SymmetricExtension{T,PL,PR,SL,SR}) where {T,PL,PR,SL,SR}
    data = unsafe_data(s)
    SymmetricExtension{T,PR,PL,SR,SL}(flipdim(data, 1), -rightindex(s))
end

"""
A compact symmetric extension is defined by a data vector `v` and corresponds to
the symmetric extension of `v` extended by zero. The symmetry can be odd or even,
and can have whole point or half point.

The compact symmetric extension is a compact sequence and it is symmetric with
respect to offset (for half point symmetry) or with respect to offset-1/2
(for whole point symmetry in which the first element of `v` is repeated).

For example, the whole point odd symmetric extension of `[1,2,3]` is the
sequence: `... 0 0 -3 -2 1 2 3 0 0 ...`
"""
struct CompactSymmetricSequence{T} <: ExtensionSequence{T}
    v           ::  Vector{T}
    odd         ::  Bool
    wholepoint  ::  Bool
    offset      ::  Int
end

CompactSymmetricSequence(v::AbstractVector, offset = firstindex(v); optional...) =
    CompactSymmetricSequence(collect(v), offset; optional...)

CompactSymmetricSequence(v::Vector, offset = 1; odd = false, wholepoint = false) =
    CompactSymmetricSequence(v, odd, wholepoint, offset)

compact_symmetric_sequence_wholepoint_even(v::AbstractVector, args...) =
    CompactSymmetricSequence(v, args...; odd=false, wholepoint=true)

compact_symmetric_sequence_halfpoint_even(v::AbstractVector, args...) =
    CompactSymmetricSequence(v, args...; odd=false, wholepoint=false)

compact_symmetric_sequence_wholepoint_odd(v::AbstractVector, args...) =
    CompactSymmetricSequence(v, args...; odd=true, wholepoint=true)

compact_symmetric_sequence_halfpoint_odd(v::AbstractVector, args...) =
    CompactSymmetricSequence(v, args...; odd=true, wholepoint=false)

unsafe_data(s::CompactSymmetricSequence) = s.v
offset(s::CompactSymmetricSequence) = s.offset

sym_odd(s::CompactSymmetricSequence) = s.odd
sym_wholepoint(s::CompactSymmetricSequence) = s.wholepoint

Base.similar(s::CompactSymmetricSequence, v::AbstractVector, offset::Int) =
    CompactSymmetricSequence(v, offset)

"Return -1 if the argument is true and 1 if it is false"
sign_parity(val::Bool) = val ? -1 : 1

mapindex(s::CompactSymmetricSequence, k) = k - offset(s) + 1
imapindex(s::CompactSymmetricSequence, l) = l + offset(s) - 1

support_rightindex(s::CompactSymmetricSequence) = offset(s) + data_length(s) - 1

support_leftindex(s::CompactSymmetricSequence) =
    offset(s) - data_length(s) + 1 + ~sym_wholepoint(s)

support(s::CompactSymmetricSequence) = support_rightindex(s), support_leftindex(s)

# We override getindex to return zero outside our embedded vector.
function Base.getindex(s::CompactSymmetricSequence, k::Int)
    v = unsafe_data(s)
    k -= offset(s)
    T = eltype(s)
    if k >= 0
        # The index is positive
        if k < length(v)
            # within the support: return the data element
            v[k+1]
        else
            # beyond the support: return zero
            zero(T)
        end
    else
        # The index is negative. We differentiate between whole and half symmetry.
        # We can simply flip the sign of the index to achieve reversal, and we
        # add one in case of whole point symmetry (so that k_flip starts at 2):
        k_flip = - k + sym_wholepoint(s)
        if k_flip > length(v)
            zero(T)
        else
            t = sign_parity(sym_odd(s))
            t*v[k_flip]
        end
    end
end

function Base.reverse(s::CompactSymmetricSequence)
    if sym_odd(s)
        if sym_wholepoint(s)
            data = -unsafe_data(s)
            data[1] = -data[1]
            CompactSymmetricSequence(data, sym_odd(s), sym_wholepoint(s), -offset(s))
        else
            CompactSymmetricSequence(-unsafe_data(s), sym_odd(s), sym_wholepoint(s), -offset(s)+1)
        end
    else
        if sym_wholepoint(s)
            CompactSymmetricSequence(copy(unsafe_data(s)), sym_odd(s), sym_wholepoint(s), -offset(s))
        else
            CompactSymmetricSequence(copy(unsafe_data(s)), sym_odd(s), sym_wholepoint(s), -offset(s)+1)
        end
    end
end

Base.ctranspose(s::CompactSymmetricSequence) = CompactSymmetricSequence(conj(reverse(s)))
