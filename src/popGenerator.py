#
# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license. See the LICENSE file in the repository root for complete license information.

from popSimulator import SimulatePopulations
import os

ONESAMP2COAL_MINALLELEFREQUENCY=0.05
mutationRate=0.012
duration_start=2
duration_range=6
missing_data_percentage=0.20
rangeNe=100,500
NeVal=200
numPOP="00256"


outputSampleSizes=(50,200)
locis=(40,320)
simulate_populations = SimulatePopulations()

for sampleSize in outputSampleSizes:
    for loci in locis:
        for i in range(1, 31):
            file_name = f"genePop{sampleSize}x{loci}_{i}"
            path = os.path.join("<<your file path>>", file_name)
            print(path)
            simulate_populations.generate_input_population(sampleSize, loci, NeVal, mutationRate, path, duration_start, duration_range, missing_data_percentage)

