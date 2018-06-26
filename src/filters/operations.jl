

"""
From a given filter h_i, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

`g_k = (-1)^k h_{pivot-k}`

The default pivot is 1.
"""
function alternating_flip(h::FIRFilter, pivot = 1)
    v = data(h)
    hflip = similar(v)
    # The first element of hflip has index pivot-rightindex(f).
    # Whether or not we need to flip its sign depend on the parity of this index:
    isodd(pivot-rightindex(h)) ? t = -1 : t = 1
    hflip[1:2:end] =  t * v[end:-2:1]
    hflip[2:2:end] = -t * v[end-1:-2:1]
    FIRFilter(hflip, pivot - rightindex(h))
end

"Compute a new 'alternating' filter satisfying `g_k = (-1)^k h_k`."
function alternating(h::FIRFilter)
    v = data(h)
    halt = similar(v)
    t = (-1)^leftindex(h)
    for k in eachindex(v)
        halt[k] = t * v[k]
        t = -t
    end
    FIRFilter(halt, leftindex(h))
end

"The first even number greater than or equal to n."
nexteven(n) = isodd(n) ? n+1 : n

"The last even number, smaller than or equal to n."
previouseven(n) = isodd(n) ? n-1 : n

"The first odd number greater than or equal to n."
nextodd(n) = isodd(n) ? n : n+1

"The previous odd number, smaller than or equal to n."
previousodd(n) = isodd(n) ? n : n-1

"Return the even part of a sequence `s`, defined by `s_e[k] = s[2k]`."
evenpart(h::FIRFilter) =
    FIRFilter([h[j] for j in nexteven(leftindex(h)):2:rightindex(h)], div(nexteven(leftindex(h)),2))

"Return the odd part of a sequence `s`, defined by `s_o[k] = s[2k+1]`."
oddpart(h::FIRFilter) =
    FIRFilter([h[j] for j in nextodd(leftindex(h)):2:rightindex(h)], div(previouseven(leftindex(h)),2))

function convolution(h1::Sequence, h2::Sequence, i)
    T = promote_type(eltype(h1), eltype(h2))
    z = zero(T)
end
