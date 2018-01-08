using MendelBase
using Search
using SearchSetup
using MendelGameteCompetition
#
# Required external modules.
#
using DataFrames    # From package DataFrames.
using Distributions # 

# @testset "initialize_optimization" begin
#     keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
#     process_keywords!(keyword, "gamete competition Control.txt", "")
#     (pedigree, person, nuclear_family, locus, snpdata,
#         locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
#         read_external_data_files(keyword)

#     loc = 1 #current locus
#     keyword["eliminate_genotypes"] = true
#     keyword["constraints"] = 1
#     keyword["goal"] = "maximize"
#     keyword["parameters"] = locus.alleles[loc]
#     # keyword["title"] = "Gamete competition analysis for " * locus.name[loc]

#     parameter = set_parameter_defaults(keyword)
#     parameter = MendelGameteCompetition.initialize_optimization_gamete_competition!(
#         locus, parameter, keyword)

#     @test parameter.travel == "search"
#     @test size(parameter.par) == (2,)
#     @test all(parameter.par .== 1.0)
#     @test size(parameter.name) == (2,)
#     @test parameter.name[1] == "tau 1   "
#     @test parameter.name[2] == "tau 2   "    
#     @test size(parameter.min) == (2,)
#     @test size(parameter.max) == (2,)
#     @test all(parameter.min .== 1e-5)
#     @test all(parameter.max .== Inf)
#     @test parameter.parameters == 2

#     @test size(parameter.constraint) == (1, 2)
#     @test parameter.constraint[1, 1] == 1.0
#     @test parameter.constraint[1, 2] == 0.0
#     @test size(parameter.constraint_level) == (1,)
#     @test parameter.constraint_level[1] == 1.0

# end

@testset "prior function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "gamete competition Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)


    loc = 1 #current locus
    keyword["eliminate_genotypes"] = true
    keyword["constraints"] = 1
    keyword["goal"] = "maximize"
    keyword["parameters"] = locus.alleles[loc]
    # keyword["title"] = "Gamete competition analysis for " * locus.name[loc]

    parameter = set_parameter_defaults(keyword)
    parameter = MendelGameteCompetition.initialize_optimization_gamete_competition!(
        locus, parameter, keyword)

    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    #
    # use dictionaries to assign prior probabilities
    #
    locus_dic = Dict()
    lucus_data = readtable("gamete competition LocusFrame.txt")
    for i in 1:length(locus.name)
        current = string(locus.name[i])

        #select the probability vector of current locus
        prob_vec = lucus_data[lucus_data[:Locus] .== current, 4] 

        #normalize the probabilities if they don't sum to 1
        if sum(prob_vec) != 1.0
            total = sum(prob_vec)
            for j in 1:length(prob_vec)
                prob_vec[j] = prob_vec[j] / total
            end
        end

        #add the probabilities to dictionary, key = order in which locus appear
        locus_dic["$i"] = prob_vec 
    end

    # loop through to compute probabilities
    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #avoids some array access errors
            if person.mother[i] != 0 continue end #prior prob doesnt exist for non founder

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                prob = MendelGameteCompetition.prior_gamete_competition(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                answer = 1.0

                # tally probabilities contributed by the 10 locus 
                for i in 1:length(locus.name)
                    prob_vec = locus_dic["$i"] #retrieve i'th locus's probability vector

                    #multi_genotype[1, i, j] is either 1 or 2, so we multiply answer by:
                    answer *= prob_vec[multi_genotype[1, i, j]] #first row ith column 
                    answer *= prob_vec[multi_genotype[2, i, j]] 

                end

                @test answer == prob
            end
        end
    end
end

@testset "transmission" begin
    
end

@testset "penetrance" begin
    
end

@testset "wrapper & basics" begin
    
end