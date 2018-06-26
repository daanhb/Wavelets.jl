print("filters ...\n")

# This is run when including test_filters.jl
function test_filters()
    test_doubly_infinite_arrays()
end

function test_doubly_infinite_arrays()
    @testset "doubly infinite sequences" begin
        @testset "zero padding" begin
            test_zero_padding()
        end
        @testset "periodic extension" begin
            test_periodic_extension()
        end
        @testset "symmetric extension" begin
            test_symmetric_extension()
        end
        @testset "compact symmetric extension" begin
            test_compact_symmetric_sequence()
        end
    end
end

function test_generic_extension_sequence_interface(z1, v)
    # Test accuracy of underlying data
    zv = [z1[i] for i in Filters.leftindex(z1):Filters.rightindex(z1)]
    vv = collect(v)
    @test sum(abs.([zv[i]-vv[i] for i in eachindex(vv)])) == 0

    # Test a shift
    z2 = shift(z1, 2)
    @test eltype(z2) == eltype(z1)
    @test Filters.leftindex(z2) == Filters.leftindex(z1)+2
    @test Filters.rightindex(z2) == Filters.rightindex(z1)+2
    @test z2[3] == z1[1]
    @test z2[4] == z1[2]

    # Test a reversal
    z1rev = reverse(z1)
    data1 = [z1[k] for k in -10:10]
    data2 = [z1rev[k] for k in 10:-1:-10]
    @test data1 == data2

    # Test conjugation
    z1conj = conj(z1)
    @test z1conj[1] == conj(z1[1])
end

function test_zero_padding()
    v1 = [1,2,3]
    # First, use a default offset and test basic functionality
    z1 = zero_padding(v1)
    @test eltype(z1) == eltype(v1)
    @test support(z1) == 1:3
    @test z1[0] == 0
    @test z1[4] == 0
    @test z1[2:3] == v1[2:3]
    test_generic_extension_sequence_interface(z1, v1)


    # Test a different offset
    offset = 2
    z3 = CompactSequence(v1, offset)
    test_generic_extension_sequence_interface(z3, v1)
    supp = Filters.joint_support(z1, z3)
    @test first(supp) == 1
    @test last(supp) == length(v1)+1

    # Try an abstract vector
    z = CompactSequence(1:5)
    @test z[3] == 3
    test_generic_extension_sequence_interface(z, 1:5)

    # Arithmetic
    z1 = CompactSequence(v1)
    z2 = CompactSequence(v1, 2)
    s1 = support(z1)
    a = z1+z2
    @test first(support(a)) == min(first(support(z1)), first(support(z2)))
    @test last(support(a)) == max(last(support(z1)), last(support(z2)))
    supp = support(a)
    @test maximum(abs.(a[supp] - z1[supp] - z2[supp])) == 0
    b = z1-z2
    @test maximum(abs.(b[supp] - z1[supp] + z2[supp])) == 0
    c = 2*z1
    @test maximum(abs.(c[s1]-2*z1[s1])) == 0
    d = -z1
    @test maximum(abs.(d[s1]+z1[s1])) == 0
    e = z1/2
    @test maximum(abs.(e[s1]-z1[s1]/2)) <= 10eps(eltype(e))
end

function test_periodic_extension()
    v1 = [1,2,3,4,5]
    z1 = PeriodicExtension(v1)
    @test eltype(z1) == eltype(v1)
    @test z1[0] == v1[end]
    @test z1[6] == v1[1]
    @test z1[2:3] == v1[2:3]
    test_generic_extension_sequence_interface(z1, v1)

    # Arithmetic
    v2 = [3,4,2,1,6]
    z1 = PeriodicExtension(v1)
    z2 = PeriodicExtension(v2, 2)
    p = period(z1)
    a = z1+z2
    @test period(a) == period(z1) == period(z2)
    s = eachindex(v2)
    @test maximum(abs.(a[s] - z1[s] - z2[s])) == 0
    b = z1-z2
    @test maximum(abs.(b[s] - z1[s] + z2[s])) == 0
    c = 2*z1
    @test maximum(abs.(c[s]-2*z1[s])) == 0
    d = -z1
    @test maximum(abs.(d[s]+z1[s])) == 0
    e = z1/2
    @test maximum(abs.(e[s]-z1[s]/2)) <= 10eps(eltype(e))
end

function test_symmetric_extension()
    v1 = [1,2,3,4,5]

    n1 = 1
    n2 = length(v1)
    z1 = symmetric_extension_wholepoint_even(v1)
    @test z1[n2+1] == z1[n2-1]
    @test z1[n2+25] == z1[n2-25]
    test_generic_extension_sequence_interface(z1, v1)


    z2 = symmetric_extension_wholepoint_odd(v1)
    @test z2[n2+1] == -z2[n2-1]
    @test z2[n2+25] == -z2[n2-25]
    test_generic_extension_sequence_interface(z2, v1)

    z3 = symmetric_extension_halfpoint_even(v1)
    @test z3[n2+1] == z3[n2]
    @test z3[n2+2] == z3[n2-1]
    @test z3[n2+25] == z3[n2-24]
    test_generic_extension_sequence_interface(z3, v1)

    z4 = symmetric_extension_halfpoint_odd(v1)
    @test z4[n2+1] == -z4[n2]
    @test z4[n2+2] == -z4[n2-1]
    @test z4[n2+25] == -z4[n2-24]
    test_generic_extension_sequence_interface(z4, v1)
end

function test_compact_symmetric_sequence()
    v1 = [1,2,3,4,5]

    n1 = 1
    n2 = length(v1)
    z1 = compact_symmetric_sequence_wholepoint_even(v1, 0)
    test_generic_extension_sequence_interface(z1, v1)
    @test z1[n2+1] == 0
    @test z1[0] == v1[1]
    @test z1[1] == z1[-1]

    z2 = compact_symmetric_sequence_wholepoint_odd(v1, 0)
    test_generic_extension_sequence_interface(z2, v1)
    @test z2[0] == v1[1]
    @test z2[n2+1] == 0
    @test z2[1] == -z2[-1]

    z3 = compact_symmetric_sequence_halfpoint_even(v1, 0)
    test_generic_extension_sequence_interface(z3, v1)
    @test z3[n2+1] == 0
    @test z3[0] == v1[1]
    @test z3[-1] == z3[0]

    z4 = compact_symmetric_sequence_halfpoint_odd(v1, 0)
    test_generic_extension_sequence_interface(z4, v1)
    @test z4[n2+1] == 0
    @test z4[0] == v1[1]
    @test z4[-1] == -z4[0]
end

test_filters()
