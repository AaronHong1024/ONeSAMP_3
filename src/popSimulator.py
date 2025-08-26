#
# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license. See the LICENSE file in the repository root for complete license information.

import msprime
import numpy as np
import random

class SimulatePopulations:

    def simulate_populations(self, sample_size, loci, neRange,  rate, duration_start, duration_range, missing_data_percentage):
        if isinstance(neRange, tuple):
            effective_population = np.random.uniform(neRange[0],neRange[1])
        else:
            effective_population = neRange
    # Simulate ancestral history and add mutations
        tree_sequence = msprime.sim_ancestry(samples=sample_size, start_time=duration_start, model=msprime.DiscreteTimeWrightFisher(duration=duration_range), ploidy=2, recombination_rate=1e-8, population_size=effective_population, sequence_length=loci, random_seed=42)
        tree_sequence = msprime.sim_mutations(tree_sequence, rate=rate, random_seed=42)

        # Convert the tree sequence to a matrix of haplotypes
        haplotypes = np.array([list(hap) for hap in tree_sequence.haplotypes()])

        # Introduce missing data
        # missing_data_proportion = 0.0
        num_sites = tree_sequence.num_sites
        num_missing = int(num_sites * missing_data_percentage * len(haplotypes)/2)

        for _ in range(num_missing):
                # Randomly choose a site and a sample
                site_index = random.randint(0, num_sites - 1)
                sample_index = random.randint(0, len(haplotypes) - 1)
                # Mark this site as missing for this sample
                haplotypes[sample_index][site_index] = 'N'  # 'N' denotes missing data

        formatted_haplotypes = []
        for i in range(0, len(haplotypes), 2):
            # Combine pairs of haplotypes into diploid individuals
            diploid = haplotypes[i].tolist() + haplotypes[i + 1].tolist()
            formatted_haplotypes.append(diploid)
        # print(effective_population)
        return formatted_haplotypes, effective_population

    # Encode haplotypes
    encoding = {'N': '00', 'A': '01', 'G': '02', 'C': '03', 'T': '04'}

    def encode_haplotypes(self, hap):
        encoded_values = [self.encoding[char] for char in hap]
        encoded_values = ' '.join([''.join(encoded_values[j:j + 2]) for j in range(0, len(encoded_values), 2)])
        return encoded_values

    # Convert simulated population to genepop format and print the effective population
    def generate_population_data(self, sample_size, loci, neRange, rate, file_name, duration_start, duration_range, missing_data_percentage):
        result = self.simulate_populations(sample_size, loci, neRange, rate, duration_start, duration_range, missing_data_percentage)
        diploid_haplotypes = result[0]
        content = "Generated genotype output\n"
        content += "\n".join(str(i+1) for i in range(loci))
        content += "\nPop"
        for i, hap in enumerate(diploid_haplotypes):
            encoded_haplotypes = self.encode_haplotypes(hap)
            content += f"\n{i+1} , {encoded_haplotypes}"
        content += f"\n{result[1]}"
        # file_name = f"genePop{sample_size}x{loci}"
        with open(file_name, 'w') as file:
            file.write(content)

    # Convert simulated population to genepop format - input population
    def generate_input_population(self, sample_size, loci, effective_population, rate, file_name, duration_start, duration_range, missing_data_percentage):
        result = self.simulate_populations(sample_size, loci, effective_population, rate, duration_start, duration_range, missing_data_percentage)
        diploid_haplotypes = result[0]
        content = "Generated genotype output\n"
        content += "\n".join(str(i+1) for i in range(loci))
        content += "\nPop"
        result = []
        for i, hap in enumerate(diploid_haplotypes):
            encoded_haplotypes = self.encode_haplotypes(hap)
            result.append(encoded_haplotypes)
            content += f"\n{i+1} , {encoded_haplotypes}"
        with open(file_name, 'w') as file:
            file.write(content)



# population = SimulatePopulations()
# population.generate_population_data(50, 20, (4,400), 0.012, "genePop50x20")
















# Print diploid haplotypes with spaces
# for i, hap in enumerate(diploid_haplotypes):
#     print(f"{i+1} , {hap}")

# tree_sequence=msprime.simulate(sample_size=50, Ne=100, length=20, mutation_rate=0.000000012, recombination_rate=3e-7)
# print(tree_sequence.genotype_matrix())

# for var in ts.variants():
#     print(var.site.position, var.alleles, var.genotypes, sep="\t")
