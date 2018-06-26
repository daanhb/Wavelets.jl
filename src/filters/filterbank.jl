# FilterBank.jl

## EXPORTS

export FilterPair, FilterBank, FilterMatrix

export polyphasematrix, modulationmatrix, polyphase_analysis!, polyphase_synthesis!,
    lowpassfilter, highpassfilter, primal_lowpassfilter, primal_highpassfilter,
    dual_lowpassfilter, dual_highpassfilter


## TYPES

"A pair of two filters."
struct FilterPair{F <: Filter}
    f1  ::  F
    f2  ::  F
end

Base.eltype(::Type{FilterPair{F}}) where F = eltype(F)

lowpassfilter(fp::FilterPair) = fp.f1
highpassfilter(fp::FilterPair) = fp.f2

function Base.getindex(fp::FilterPair, k::Int)
    if k == 1
        fp.f1
    elseif k == 2
        fp.f2
    else
        throw(BoundsError())
    end
end

Base.ctranspose(fb::FilterPair) = FilterPair(fb.f1', fb.f2')


"A 2x2 matrix of filters."
struct FilterMatrix{F <: Filter}
    a11 ::  F
    a12 ::  F
    a21 ::  F
    a22 ::  F
end

Base.eltype(::Type{FilterMatrix{F}}) where F = eltype(F)

Base.transpose(m::FilterMatrix) = FilterMatrix(m.a11, m.a21, m.a12, m.a22)

Base.ctranspose(m::FilterMatrix) = FilterMatrix(m.a11', m.a21', m.a12', m.a22')

(m::FilterMatrix)(z) = [ztransform(m.a11, z) ztransform(m.a12, z); ztransform(m.a21, z) ztransform(m.a22, z)]

polyphasematrix(fp::FilterPair) = FilterMatrix(evenpart(fp[1]), evenpart(fp[2]), oddpart(fp[1]), oddpart(fp[2]))

modulationmatrix(fp::FilterPair) = FilterMatrix(fp[1], alternating(fp[1]), fp[2], alternating(fp[2]))



"A FilterBank groups several objects related to a two-phase FilterBank."
struct FilterBank{F <: Filter}
    "The primal filter pair, to be used on the synthesis side."
    primal_pair     ::  FilterPair{F}
    "The dual filter pair, to be used on the analysis side."
    dual_pair       ::  FilterPair{F}
    "The polyphase matrix for the synthesis side."
    pm_synthesis    ::  FilterMatrix{F}
    "The polyphase matrix for the analysis side."
    pm_analysis     ::  FilterMatrix{F}
end

FilterBank(lowpass::Filter) = FilterBank(FilterPair(lowpass, alternating_flip(lowpass)))

FilterBank(primal_lowpass::Filter, dual_lowpass::Filter) =
    FilterBank(FilterPair(primal_lowpass, alternating_flip(dual_lowpass)),
        FilterPair(dual_lowpass, alternating_flip(primal_lowpass)))

# If no dual pair is given, assume orthogonality.
# TODO: verify orthogonality
FilterBank(primal_pair::FilterPair) = FilterBank(primal_pair, primal_pair)

FilterBank(primal_pair::FilterPair, dual_pair::FilterPair) =
    FilterBank(primal_pair, dual_pair, polyphasematrix(primal_pair), polyphasematrix(dual_pair)')

Base.eltype(::Type{FilterBank{F}}) where F = eltype(F)

primal_lowpassfilter(fb::FilterBank) = lowpassfilter(fb.primal_pair)
primal_highpassfilter(fb::FilterBank) = highpassfilter(fb.primal_pair)

dual_lowpassfilter(fb::FilterBank) = lowpassfilter(fb.dual_pair)
dual_highpassfilter(fb::FilterBank) = highpassfilter(fb.dual_pair)

# Split the (finite) signal x embedded in (infinite) embedding into a low pass signal y1
# and a high pass signal y2. The filters are given by the polyphasematrix f.
function polyphase_analysis!(y1, y2, x, f::FilterMatrix, embedding)
    Heven = f.a11
    Hodd = f.a12
    Geven = f.a21
    Godd = f.a22

    T = eltype(y1)
    # Determine for which range the convolution of the filters with x results
    # in evaluations of x that not occur outside of the embedding.
    # Remark (TODO?): this could be even more efficient.
    L = length(x)
    lower_i = min(length(y1)-1, maximum(map(lastindex,(Heven,Geven,Godd,Hodd))))
    L%2 == 0?
      upper_i = L>>1 - 1 + minimum(map(firstindex,(Heven,Geven,Godd,Hodd))) :
      upper_i = L>>1 - L%2 + minimum(map(firstindex,(Heven,Geven,Godd,Hodd)))
    upper_i = max(0, upper_i)
    # Lower boundary region
    @inbounds for i in 0:lower_i
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * embedding[x, 2*i-2*l]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * embedding[x, 2*i-2*l]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * embedding[x, 2*i-2*l+1]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * embedding[x, 2*i-2*l+1]
        end
        y1[i+1] = y1i
        y2[i+1] = y2i
    end

    # Middle region
    @inbounds for i in lower_i+1:upper_i-1
      # for l in minimum(map(firstindex,(Heven,Geven,Hodd,Godd))):maximum(map(lastindex,(Heven,Geven,Hodd,Godd)))
      #   @assert (0<=2(i-l)<L) && (0<=2(i-l)+1<L)
      # end
      # @assert 0<=i<min(length(y1),length(y2))
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * x[2*i-2*l+1]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * x[2*i-2*l+1]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * x[2*i-2*l+2]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * x[2*i-2*l+2]
        end
        y1[i+1] = y1i
        y2[i+1] = y2i
    end

    # Upper boundary region
    @inbounds for i in upper_i:length(y1)-1
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * embedding[x, 2*i-2*l]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * embedding[x, 2*i-2*l]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * embedding[x, 2*i-2*l+1]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * embedding[x, 2*i-2*l+1]
        end
        y1[i+1] = y1i
        if i+1 <= length(y2)
            y2[i+1] = y2i
        end
    end
end

function polyphase_synthesis!(x, y1, y2, f::FilterMatrix, embedding)
    Heven = f.a11
    Geven = f.a12
    Hodd = f.a21
    Godd = f.a22

    T = eltype(x)
    l1 = length(y1)
    l2 = length(y2)

    for j in 0:length(x)>>1 - 1
        xj_e = zero(T)
        xj_o = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            xj_e += Heven[l] * embedding[y1,j-l, l1]
        end
        for l = firstindex(Geven):lastindex(Geven)
            xj_e += Geven[l] * embedding[y2,j-l, l2]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            xj_o += Hodd[l] * embedding[y1,j-l, l1]
        end
        for l = firstindex(Godd):lastindex(Godd)
            xj_o += Godd[l] * embedding[y2,j-l, l2]
        end
        x[2*j+1] = xj_e
        x[2*j+2] = xj_o
    end
end

function polyphase_analysis!(y::AbstractVector{T}, l1::Int, l2::Int, x::AbstractVector{T}, L::Int, f::FilterMatrix, embedding) where T
    Heven = f.a11
    Hodd = f.a12
    Geven = f.a21
    Godd = f.a22

    # Determine for which range the convolution of the filters with x results
    # in evaluations of x that not occur outside of the embedding.
    # Remark (TODO?): this could be even more efficient.


    lower_i = min(l1-1, maximum(map(lastindex,(Heven,Geven,Godd,Hodd))))
    L%2 == 0?
      upper_i = L>>1 - 1 + minimum(map(firstindex,(Heven,Geven,Godd,Hodd))) :
      upper_i = L>>1 - L%2 + minimum(map(firstindex,(Heven,Geven,Godd,Hodd)))
    upper_i = max(0, upper_i)

    # Lower boundary region
    @inbounds for i in 0:lower_i
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * embedding[x, 2*i-2*l, L]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * embedding[x, 2*i-2*l, L]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * embedding[x, 2*i-2*l+1, L]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * embedding[x, 2*i-2*l+1, L]
        end
        y[i+1] = y1i
        y[l1+i+1] = y2i
    end

    # Middle region
    @inbounds for i in lower_i+1:upper_i-1
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * x[2*i-2*l+1]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * x[2*i-2*l+1]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * x[2*i-2*l+2]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * x[2*i-2*l+2]
        end
        y[i+1] = y1i
        y[l1+i+1] = y2i
    end

    # Upper boundary region
    @inbounds for i in upper_i:l1-1
        y1i = zero(T)
        y2i = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            y1i += Heven[l] * embedding[x, 2*i-2*l, L]
        end
        for l = firstindex(Geven):lastindex(Geven)
            y2i += Geven[l] * embedding[x, 2*i-2*l, L]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            y1i += Hodd[l] * embedding[x, 2*i-2*l+1, L]
        end
        for l = firstindex(Godd):lastindex(Godd)
            y2i += Godd[l] * embedding[x, 2*i-2*l+1, L]
        end
        y[i+1] = y1i
        if i+1 <= l2
            y[l1+i+1] = y2i
        end
    end
end

function polyphase_synthesis!(x, L::Int, y, l1::Int, l2::Int, f::FilterMatrix, embedding)
    Heven = f.a11
    Geven = f.a12
    Hodd = f.a21
    Godd = f.a22

    T = eltype(x)
    l = L>>1 - 1

    # TODO: distinguish between boundary regions and the middle region
    @inbounds for j in 0:l
        xj_e = zero(T)
        xj_o = zero(T)
        for l = firstindex(Heven):lastindex(Heven)
            xj_e += Heven[l] * embedding[y,j-l, l1]
        end
        for l = firstindex(Geven):lastindex(Geven)
            xj_e += Geven[l] * embedding[y,j-l, l2, l1]
        end
        for l = firstindex(Hodd):lastindex(Hodd)
            xj_o += Hodd[l] * embedding[y,j-l, l1]
        end
        for l = firstindex(Godd):lastindex(Godd)
            xj_o += Godd[l] * embedding[y,j-l, l2, l1]
        end
        x[(j<<1)+1] = xj_e
        x[(j<<1)+2] = xj_o
    end
end
