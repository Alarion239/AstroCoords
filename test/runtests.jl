using AstroCoords
using Test
using Random

@testset "AstroCoords.jl" begin
    @testset "Constructor Tests" begin
        @testset "Basic Construction" begin
            # Basic Spherical
            s = Spherical(π/4, π/6)
            @test s.latitude ≈ π/4
            @test s.longitude ≈ π/6
            @test lat(s) ≈ π/4
            @test lon(s) ≈ π/6
            
            # Basic SphericalD
            sd = SphericalD(π/4, π/6, 2.0)
            @test sd.latitude ≈ π/4
            @test sd.longitude ≈ π/6
            @test sd.distance ≈ 2.0
            @test lat(sd) ≈ π/4
            @test lon(sd) ≈ π/6
            @test dist(sd) ≈ 2.0
            
            # Basic Cartesian
            c = Cartesian(1.0, 0.0, 0.0)
            @test norm([c.x, c.y, c.z]) ≈ 1.0  # Should be normalized
            
            # Basic CartesianD
            cd = CartesianD(1.0, 2.0, 3.0)
            @test cd.x ≈ 1.0
            @test cd.y ≈ 2.0
            @test cd.z ≈ 3.0
            @test dist(cd) ≈ sqrt(14)
        end
        
        @testset "Type Promotion" begin
            # Mixed numeric types
            s = Spherical(1, 2.0)  # Int and Float64
            @test typeof(s.latitude) == typeof(s.longitude) == Float64
            
            sd = SphericalD(1, 2.0, 3.0f0)  # Mixed types for angles, separate for distance
            @test typeof(sd.latitude) == typeof(sd.longitude) == Float64
            @test typeof(sd.distance) == Float32
            
            c = Cartesian(1, 2.0, 3.0f0)
            @test eltype(typeof(c)) == Float64  # Should promote to highest precision
        end
        
        @testset "Constructor Edge Cases" begin
            # Zero vector should throw error for Cartesian
            @test_throws ArgumentError Cartesian(0.0, 0.0, 0.0)
            
            # Very small vector should still work
            tiny = 1e-15
            c = Cartesian(tiny, 0.0, 0.0)
            @test norm(c) ≈ 1.0
            
            # Large coordinates
            large = 1e10
            cd = CartesianD(large, large, large)
            @test cd.x == large
            
            # Poles (mathematical edge case)
            pole_north = Spherical(π/2, 0.0)
            pole_south = Spherical(-π/2, 0.0)
            @test lat(pole_north) ≈ π/2
            @test lat(pole_south) ≈ -π/2
        end
    end
    
    @testset "Conversion Tests" begin
        @testset "Basic Conversions" begin
            # Spherical to Cartesian
            s = Spherical(0.0, 0.0)  # Equator, prime meridian
            c = Cartesian(s)
            @test c.x ≈ 1.0 atol=1e-15
            @test c.y ≈ 0.0 atol=1e-15
            @test c.z ≈ 0.0 atol=1e-15
            
            # Cartesian to Spherical
            c2 = Cartesian(1.0, 0.0, 0.0)
            s2 = Spherical(c2)
            @test lat(s2) ≈ 0.0 atol=1e-15
            @test lon(s2) ≈ 0.0 atol=1e-15
            
            # Round-trip conversion
            original = Spherical(π/4, π/3)
            round_trip = Spherical(Cartesian(original))
            @test isapprox(original, round_trip, atol=1e-15)
        end
        
        @testset "Distance Conversions" begin
            # SphericalD to CartesianD
            sd = SphericalD(π/4, π/6, 2.0)
            cd = CartesianD(sd)
            @test dist(cd) ≈ 2.0
            
            # CartesianD to SphericalD
            cd2 = CartesianD(1.0, 1.0, 1.0)
            sd2 = SphericalD(cd2)
            @test dist(sd2) ≈ sqrt(3)
            
            # Round-trip with distance
            original_d = SphericalD(π/3, π/4, 5.0)
            round_trip_d = SphericalD(CartesianD(original_d))
            @test isapprox(original_d, round_trip_d, atol=1e-14)
        end
        
        @testset "Conversion Edge Cases" begin
            # North pole
            north_pole = Spherical(π/2, 0.0)
            c_north = Cartesian(north_pole)
            @test c_north.x ≈ 0.0 atol=1e-15
            @test c_north.y ≈ 0.0 atol=1e-15
            @test c_north.z ≈ 1.0 atol=1e-15
            
            # South pole
            south_pole = Spherical(-π/2, 0.0)
            c_south = Cartesian(south_pole)
            @test c_south.z ≈ -1.0 atol=1e-15
            
            # Longitude wraparound (should be consistent)
            s1 = Spherical(π/4, 0.0)
            s2 = Spherical(π/4, 2π)
            c1 = Cartesian(s1)
            c2 = Cartesian(s2)
            @test isapprox(c1.x, c2.x, atol=1e-15)
            @test isapprox(c1.y, c2.y, atol=1e-15)
        end
    end
    
    @testset "Equality Tests" begin
        @testset "Same Type Equality" begin
            s1 = Spherical(π/4, π/6)
            s2 = Spherical(π/4, π/6)
            s3 = Spherical(π/4, π/5)
            
            @test s1 == s2
            @test s1 != s3
            @test isequal(s1, s2)
            @test !isequal(s1, s3)
        end
        
        @testset "Cross-Type Equality" begin
            s = Spherical(π/4, π/6)
            c = Cartesian(s)
            
            @test s == c
            @test c == s
            @test isequal(s, c)
            @test isequal(c, s)
        end
        
        @testset "Distance vs Unit Sphere Equality" begin
            s = Spherical(π/4, π/6)
            sd_unit = SphericalD(π/4, π/6, 1.0)
            sd_non_unit = SphericalD(π/4, π/6, 2.0)
            
            @test s == sd_unit  # Unit distance should equal unit sphere
            @test s != sd_non_unit  # Non-unit distance should not equal unit sphere
            @test isequal(s, sd_unit)
            @test !isequal(s, sd_non_unit)
        end
        
        @testset "NaN Handling" begin
            s_nan1 = Spherical(NaN, π/6)
            s_nan2 = Spherical(NaN, π/6)
            s_normal = Spherical(π/4, π/6)
            
            # NaN should not equal itself with ==
            @test !(s_nan1 == s_nan2)
            @test !(s_nan1 == s_normal)
            
            # But should with isequal
            @test isequal(s_nan1, s_nan2)
            @test !isequal(s_nan1, s_normal)
        end
    end
    
    @testset "Hash Consistency Tests" begin
        @testset "Basic Hash Consistency" begin
            s = Spherical(π/4, π/6)
            c = Cartesian(s)
            
            # Equal objects must have equal hashes
            @test hash(s) == hash(c)
            
            # Same representation should have same hash
            s2 = Spherical(π/4, π/6)
            @test hash(s) == hash(s2)
        end
        
        @testset "Set Operations" begin
            s = Spherical(π/4, π/6)
            c = Cartesian(s)
            
            # Should work correctly in Sets
            set = Set([s, c])
            @test length(set) == 1  # Should be treated as same element
            
            # Dict operations
            dict = Dict(s => "test")
            @test dict[c] == "test"  # Should find same key
        end
        
        @testset "Distance Hash Distinction" begin
            s = Spherical(π/4, π/6)
            sd = SphericalD(π/4, π/6, 2.0)
            
            # Different distance representations should have different hashes
            # (even if they represent same direction)
            @test hash(s) != hash(sd)
        end
    end
    
    @testset "Approximate Equality Tests" begin
        @testset "Basic Approximate Equality" begin
            s1 = Spherical(π/4, π/6)
            s2 = Spherical(π/4 + 1e-15, π/6)
            
            @test !( s1 == s2)  # Should not be exactly equal
            @test isapprox(s1, s2)  # But should be approximately equal
        end
        
        @testset "Cross-Type Approximate Equality" begin
            s = Spherical(π/4, π/6)
            c = Cartesian(s)
            
            # Add tiny perturbation
            c_perturbed = Cartesian(c.x + 1e-15, c.y, c.z)  # This will get normalized
            
            @test isapprox(s, c_perturbed, atol=1e-10)
        end
        
        @testset "Custom Tolerances" begin
            s1 = Spherical(1.0, 2.0)
            s2 = Spherical(1.1, 2.0)
            
            @test !isapprox(s1, s2)  # Default tolerance
            @test isapprox(s1, s2, atol=0.2)  # Large absolute tolerance
            @test isapprox(s1, s2, rtol=0.1)  # Large relative tolerance
        end
        
        @testset "NaN in Approximate Equality" begin
            s_nan = Spherical(NaN, π/6)
            s_normal = Spherical(π/4, π/6)
            
            @test !isapprox(s_nan, s_normal)
            @test isapprox(s_nan, s_nan, nans=true)  # NaN equal to itself with nans=true
        end
    end
    
    @testset "Utility Functions Tests" begin
        @testset "Distance Function" begin
            # Distance between identical points
            s1 = Spherical(π/4, π/6)
            @test distance(s1, s1) ≈ 0.0 atol=1e-15
            
            # Distance between opposite points on sphere
            north = Spherical(π/2, 0.0)
            south = Spherical(-π/2, 0.0)
            @test distance(north, south) ≈ π atol=1e-15
            
            # Distance between points 90 degrees apart
            p1 = Spherical(0.0, 0.0)
            p2 = Spherical(π/2, 0.0)
            @test distance(p1, p2) ≈ π/2 atol=1e-15
        end
        
        @testset "Norm Function" begin
            # Unit vector
            c = Cartesian(1.0, 0.0, 0.0)
            @test norm(c) ≈ 1.0
            
            # Non-unit vector
            cd = CartesianD(3.0, 4.0, 0.0)
            @test norm(cd) ≈ 5.0
        end
    end
    
    @testset "Edge Cases and Special Values" begin
        @testset "Infinity Handling" begin
            # Infinite coordinates in CartesianD
            cd_inf = CartesianD(Inf, 1.0, 1.0)
            @test isinf(dist(cd_inf))
            
            # Should not break conversions (though result may be inf/nan)
            sd_from_inf = SphericalD(cd_inf)
            @test isinf(dist(sd_from_inf))
        end
        
        @testset "Very Small/Large Values" begin
            # Very small coordinates
            tiny = 1e-100
            s_tiny = Spherical(tiny, tiny)
            c_tiny = Cartesian(s_tiny)
            @test norm(c_tiny) ≈ 1.0
            
            # Very large coordinates in SphericalD
            huge_dist = 1e100
            sd_huge = SphericalD(π/4, π/6, huge_dist)
            cd_huge = CartesianD(sd_huge)
            @test dist(cd_huge) ≈ huge_dist
        end
        
        @testset "Precision at Poles" begin
            # North pole with different longitudes should give same Cartesian result
            north1 = Spherical(π/2, 0.0)
            north2 = Spherical(π/2, π)
            c1 = Cartesian(north1)
            c2 = Cartesian(north2)
            
            @test isapprox(c1.x, c2.x, atol=1e-15)
            @test isapprox(c1.y, c2.y, atol=1e-15)
            @test isapprox(c1.z, c2.z, atol=1e-15)
        end
    end
    
    @testset "Random Property Tests" begin
        Random.seed!(42)  # For reproducible tests
        
        @testset "Conversion Round-Trip Properties" begin
            for _ in 1:100
                # Random spherical coordinates
                lat = (rand() - 0.5) * π  # -π/2 to π/2
                lon = rand() * 2π  # 0 to 2π
                dist_val = rand() * 100 + 0.1  # Avoid zero distance
                
                # Test Spherical round-trip
                s_orig = Spherical(lat, lon)
                s_round = Spherical(Cartesian(s_orig))
                @test isapprox(s_orig, s_round, atol=1e-14)
                
                # Test SphericalD round-trip
                sd_orig = SphericalD(lat, lon, dist_val)
                sd_round = SphericalD(CartesianD(sd_orig))
                @test isapprox(sd_orig, sd_round, atol=1e-13)
            end
        end
        
        @testset "Distance Preservation" begin
            for _ in 1:50
                lat = (rand() - 0.5) * π
                lon = rand() * 2π
                dist_val = rand() * 10 + 0.1
                
                sd = SphericalD(lat, lon, dist_val)
                cd = CartesianD(sd)
                
                # Distance should be preserved in conversion
                @test isapprox(dist(sd), dist(cd), rtol=1e-14)
            end
        end
        
        @testset "Equality Consistency" begin
            for _ in 1:50
                lat = (rand() - 0.5) * π
                lon = rand() * 2π
                
                s = Spherical(lat, lon)
                c = Cartesian(s)
                
                # Cross-type equality should be consistent
                @test s == c
                @test c == s
                @test isequal(s, c)
                @test hash(s) == hash(c)
            end
        end
    end
    
    @testset "Pretty Printing Tests" begin
        # Test that show methods don't crash
        s = Spherical(π/4, π/6)
        sd = SphericalD(π/4, π/6, 2.0)
        c = Cartesian(1.0, 0.0, 0.0)
        cd = CartesianD(1.0, 2.0, 3.0)
        
        @test sprint(show, s) isa String
        @test sprint(show, sd) isa String
        @test sprint(show, c) isa String
        @test sprint(show, cd) isa String
        
        # Test that type information is shown for SphericalD
        sd_mixed = SphericalD(1.0f0, 2.0f0, 3.0)  # Float32 angles, Float64 distance
        output = sprint(show, sd_mixed)
        @test occursin("Float32", output)
        @test occursin("Float64", output)
    end
    
end
