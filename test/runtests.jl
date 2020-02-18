include("common.jl")

macro include_testset(filename)
    @assert filename isa AbstractString
    quote
        @testset $(filename) begin
            include($(filename))
        end
    end
end

####
#### unit tests
####

@include_testset("test-theta.jl")
@include_testset("test-weierstrass.jl")
