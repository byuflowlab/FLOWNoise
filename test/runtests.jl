import FLOWNoise
using Base.Test
include("testingfunctions.jl")

@testset "Post-processing: Time to frequency domain" begin
    @testset "Case 1" begin # currently only converts time-domain output to frequency-domain output and compares to WOPWOP frequency domain

        filepath = "./data/case1/"
        windowtype = "None"

        thickness_dB_calc, loading_dB_calc, total_dB_calc,
            thickness_dBA_calc, loading_dBA_calc, total_dBA_calc,
            thickness_dB_actual, loading_dB_actual, total_dB_actual,
            thickness_dBA_actual, loading_dBA_actual, total_dBA_actual,
            thickness_OASPL_calc, loading_OASPL_calc, total_OASPL_calc, 
            thickness_OASPLA_calc, loading_OASPLA_calc, total_OASPLA_calc, 
            thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual, 
            thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual = prepare_test(filepath, windowtype)

        tol_case1 = 1e-2 # 1% error

        @testset "OASPL" begin
            @test isapprox(percent_error(thickness_OASPL_actual[1], thickness_OASPL_calc), 0.0, atol=tol_case1) 
            @test isapprox(percent_error(loading_OASPL_actual[1], loading_OASPL_calc), 0.0, atol=tol_case1) 
            @test isapprox(percent_error(total_OASPL_actual[1], total_OASPL_calc), 0.0, atol=tol_case1) 
        end 
        @testset "OASPLA" begin
            @test isapprox(percent_error(thickness_OASPLA_actual[1], thickness_OASPLA_calc), 0.0, atol=tol_case1) 
            @test isapprox(percent_error(loading_OASPLA_actual[1], loading_OASPLA_calc), 0.0, atol=tol_case1) 
            @test isapprox(percent_error(total_OASPLA_actual[1], total_OASPLA_calc), 0.0, atol=tol_case1)   
        end
        @testset "Thickness SPL" begin #for loops are split up to make the tests more distinct, easier to see what's wrong instead of lumping together
            for i in 1:length(thickness_dB_actual)
                @test isapprox(percent_error(thickness_dB_actual[i], thickness_dB_calc[i]), 0.0, atol=tol_case1)
            end
        end
        @testset "Loading SPL" begin
            for i in 1:length(thickness_dB_actual)
                @test isapprox(percent_error(loading_dB_actual[i], loading_dB_calc[i]), 0.0, atol=tol_case1)
            end
        end
        @testset "Total SPL" begin
            for i in 1:length(thickness_dB_actual)
                @test isapprox(percent_error(total_dB_actual[i], total_dB_calc[i]), 0.0, atol=tol_case1)
            end
        end
        @testset "Thickness SPLA" begin
            for i in 1:length(thickness_dBA_actual)
                @test isapprox(percent_error(thickness_dBA_actual[i], thickness_dBA_calc[i]), 0.0, atol=tol_case1)
            end
        end
        @testset "Loading SPLA" begin
            for i in 1:length(thickness_dBA_actual)
                @test isapprox(percent_error(loading_dBA_actual[i], loading_dBA_calc[i]), 0.0, atol=tol_case1)
            end
        end
        @testset "Total SPLA" begin
            for i in 1:length(thickness_dBA_actual)
                @test isapprox(percent_error(total_dBA_actual[i], total_dBA_calc[i]), 0.0, atol=tol_case1)
            end
        end
    end

    # @testset "Case 2" begin #- Look at low/high pass filtering

    #     filepath = "./data/case2/"
    #     windowtype = "Hann"

    #     thickness_OASPL_calc, loading_OASPL_calc, total_OASPL_calc, 
    #         thickness_OASPLA_calc, loading_OASPLA_calc, total_OASPLA_calc, 
    #         thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual, 
    #         thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual = prepare_test(filepath, windowtype)

    #     tol_case1 = 1e-2

    #     @test isapprox(OASPL_thickness_dB, thickness_OASPL[1], atol=tol_case1) 
    #     @test isapprox(OASPL_loading_dB, loading_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dB, total_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_thickness_dBA, thickness_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_loading_dBA, loading_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dBA, total_OASPLA[1], atol=tol_case1)


    # end

    # @testset "Case 4" begin - need to look at filtering

    # filepath = "./data/case5/"
    #     windowtype = "None"

    #     thickness_OASPL_calc, loading_OASPL_calc, total_OASPL_calc, 
    #         thickness_OASPLA_calc, loading_OASPLA_calc, total_OASPLA_calc, 
    #         thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual, 
    #         thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual = prepare_test(filepath, windowtype)

    #     tol_case1 = 1e-2

    #     @test isapprox(OASPL_thickness_dB, thickness_OASPL[1], atol=tol_case1) 
    #     @test isapprox(OASPL_loading_dB, loading_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dB, total_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_thickness_dBA, thickness_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_loading_dBA, loading_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dBA, total_OASPLA[1], atol=tol_case1)


    # end

    # @testset "Case 5" begin - SPL found only for two specific ranges, not each frequency in the spectrum

    # filepath = "./data/case5/"
    #     windowtype = "None"

    #     thickness_OASPL_calc, loading_OASPL_calc, total_OASPL_calc, 
    #         thickness_OASPLA_calc, loading_OASPLA_calc, total_OASPLA_calc, 
    #         thickness_OASPL_actual, loading_OASPL_actual, total_OASPL_actual, 
    #         thickness_OASPLA_actual, loading_OASPLA_actual, total_OASPLA_actual = prepare_test(filepath, windowtype)

    #     tol_case1 = 1e-2

    #     @test isapprox(OASPL_thickness_dB, thickness_OASPL[1], atol=tol_case1) 
    #     @test isapprox(OASPL_loading_dB, loading_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dB, total_OASPL[1], atol=tol_case1)
    #     @test isapprox(OASPL_thickness_dBA, thickness_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_loading_dBA, loading_OASPLA[1], atol=tol_case1)
    #     @test isapprox(OASPL_total_dBA, total_OASPLA[1], atol=tol_case1)

    # end

end