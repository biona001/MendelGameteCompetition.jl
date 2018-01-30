using MendelBase
using Search
using SearchSetup
using MendelGameteCompetition
#
# Required external modules.
#
using DataFrames    # From package DataFrames.
using Distributions # 

@testset "initialize_optimization" begin
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

    @test parameter.travel == "search"
    @test size(parameter.par) == (2,)
    @test all(parameter.par .== 1.0)
    @test size(parameter.name) == (2,)
    @test parameter.name[1] == "tau 1   " #not sure why there's all these extra spaces
    @test parameter.name[2] == "tau 2   "    
    @test size(parameter.min) == (2,)
    @test size(parameter.max) == (2,)
    @test all(parameter.min .== 1e-5)
    @test all(parameter.max .== Inf)
    @test parameter.parameters == 2

    @test size(parameter.constraint) == (1, 2)
    @test parameter.constraint[1, 1] == 1.0
    @test parameter.constraint[1, 2] == 0.0
    @test size(parameter.constraint_level) == (1,)
    @test parameter.constraint_level[1] == 1.0

end

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
    keyword["title"] = "Gamete competition analysis for " * locus.name[loc]

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

                    #multi_genotype[1, i, j] is either 0 or 1 or 2. I'm not sure why
                    # it can be 0, but if it is, skip it. 
                    # Otherwise we multiply answer by:
                    if multi_genotype[1, i, j] == 0 || multi_genotype[2, i, j] == 0
                        continue
                    end
                    answer *= prob_vec[multi_genotype[1, i, j]] #first row ith column 
                    answer *= prob_vec[multi_genotype[2, i, j]] 
                end

                @test answer == prob
            end
        end
    end
end

@testset "transmission" begin
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
    keyword["title"] = "Gamete competition analysis for " * locus.name[loc]

    parameter = set_parameter_defaults(keyword)
    parameter = MendelGameteCompetition.initialize_optimization_gamete_competition!(
        locus, parameter, keyword)

    par = parameter.par #this is [1.0; 1.0] under the null hypothesis of mendelian segregation 
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            if operation != transmission_array; continue; end #this avoids some array access errors

            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            j = instruction.extra[n][4]

            # need 2 genotypes to run transmission, so we construct them
            # and give them the same names as used in elston_stewart_evaluation
            i_genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                i_genotypes, i)

            maternal = !person.male[i]
            j_genotypes = MendelBase.genotype_count(person, locus, j, start, finish)
            gamete = MendelBase.construct_gametes(person, locus, start, finish, j_genotypes, j,
                                         maternal)

            # loop through all possible genotypes, and see if transmission probability is correct
            for l = 1:j_genotypes
                for k = 1:i_genotypes
                    trans = MendelGameteCompetition.transmission_gamete_competition(person, 
                        locus, gamete[:, l], multi_genotype[:, :, k], par, keyword, 
                        start, finish, i, j)

                    #
                    #compute the probaiblity that a multi_genotype pass down the gamete
                    #

                    #find first instance where locus allele is not 0
                    first_nonzero_locus_gamete = findfirst(gamete[:, l])
                    first_nonzero_locus_mom = findfirst(multi_genotype[1, :, k])
                    first_nonzero_locus_dad = findfirst(multi_genotype[2, :, k])

                    @test first_nonzero_locus_gamete == first_nonzero_locus_mom == 
                        first_nonzero_locus_dad

                    offspring_gamete = gamete[first_nonzero_locus_gamete, l]
                    mom = multi_genotype[1, first_nonzero_locus_mom, k] #first non-zero entry 
                    dad = multi_genotype[2, first_nonzero_locus_dad, k]

                    if offspring_gamete == mom && offspring_gamete == dad
                        @test trans == 1.0
                    elseif offspring_gamete == mom || offspring_gamete == dad
                        # if parent heterozygous, then transmit allele with probability 
                        # 1 / (1 + 1) = 0.5
                        # where 1 came from par = [1.0; 1.0] because Mendelian segregation 
                        # corresponds to the choice of τ_i = 1 for all i
                        @test trans == 0.5
                    else 
                        @test trans == 0.0
                    end

                end
            end
        end    
    end
end

@testset "penetrance function" begin
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

    parameter = set_parameter_defaults(keyword)

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)
    par = parameter.par #this is [0.0]

    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #this avoids some array access errors

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                number = MendelGameteCompetition.penetrance_gamete_competition(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                @test number == 1.0 #this is always 1 
            end
        end
    end
end

@testset "are wrappers doing what they should be doing" begin
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
    keyword["gamete_competition_table"] = "gamete competition Table Output.txt"

    skipped_loci = MendelGameteCompetition.gamete_competition_option(pedigree, person, 
        nuclear_family, locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
    result = GameteCompetition("gamete competition Control.txt")

    @test skipped_loci == 0
    @test result == nothing #returning nothing implies no error have been thrown
end

#current coverage = (95, 105) ≈ 90.5%







