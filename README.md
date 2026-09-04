[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A Nested Branch-and-Cut Approach for Robust Freight Service Planning on Urban Rail Transit Networks

This archive is prepared for distribution in association with the [INFORMS
Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the
[MIT License](LICENSE.txt).

The software and data in this repository are a snapshot of the software and
data used in the research reported in *A Nested Branch-and-Cut Approach for 
Robust Freight Service Planning on Urban Rail Transit Networks* by Tao Feng, 
Qinghe Sun, Shuaian Wang, and Lingxiao Wu.

## Cite

The article DOI and the archival repository DOI will be added after formal
publication. Until then, please cite the associated manuscript and this
repository using the following BibTeX record.

```bibtex
@misc{Feng2026ROFSP,
  author    = {Tao Feng and Qinghe Sun and Shuaian Wang and Lingxiao Wu},
  publisher = {INFORMS Journal on Computing},
  title     = {A Nested Branch-and-Cut
               Approach for Robust Freight Service Planning on Urban Rail
               Transit Networks},
  year      = {2026},
  url       = {https://github.com/fengtao818/2025.1492.07},
  note      = {Computational replication package},
}
```

## Description

- This repository is the computational replication package for the associated
  manuscript.
- It provides the C\# implementation of the nested branch-and-cut (NeBAC)
  approach and the comparison configurations used in the computational study.
- It includes the curated instance data, disruption scenarios, and the final
  workbook containing the numerical results reported in the manuscript and its
  Electronic Companion.
- For questions about the package, please contact the corresponding author
  listed in [AUTHORS.txt](AUTHORS.txt).

## Repository Structure

The repository is organized as follows.

### DATASETS

- The instance data and disruption scenarios are provided in
  [`Results/UploadData`](Results/UploadData).
- Each structural instance is organized as `N**D**/data1` through `data5`.
- The independent (Bernoulli-based) and line-based 1,000-scenario files are
  stored in `1000_instances`.
- For a detailed explanation of the instance layout and the mapping between
  inputs and manuscript tables, see
  [`Results/UploadData/README.md`](Results/UploadData/README.md).

### RESULTS

- The final computational results are provided in
  [`Results/UploadData/resutls/RobustTrainresults.xlsx`](Results/UploadData/resutls/RobustTrainresults.xlsx).
- The workbook contains the results reported in the main text and the
  Electronic Companion, organized by manuscript section.
- See [`Results/README.md`](Results/README.md) for the interpretation of each
  worksheet and its relation to the reported tables.

### CODE

- The implementation is provided in
  [`ROernightTrainLogistics`](ROernightTrainLogistics).
- `RoOvernightTrainLogistics.sln` is the Visual Studio solution;
  `RoOvernightTrainLogistics/Program.cs` is the program entry point; and
  `RoOvernightTrainLogistics/ColBDGen.cs` implements the main decomposition
  and cut-generation procedures.
- The implementation was developed using Visual Studio 2019, .NET Framework
  4.8, and IBM ILOG CPLEX Optimization Studio 22.1.1. CPLEX is proprietary and
  is not distributed with this repository.
- The code reads `Input_Nodes.csv`, `Input_Links.csv`, `Input_Paths.csv`, and
  `Input_Disruptions.csv` from `ROernightTrainLogistics/Dataset/` and writes
  run-specific output to the repository-level `TestLog/` directory.
- Detailed installation, dependency, and execution instructions are available
  in [`ROernightTrainLogistics/README.md`](ROernightTrainLogistics/README.md).

## Replicating

To run the NeBAC implementation for one instance:

1. Install Visual Studio 2019, .NET Framework 4.8, and a licensed copy of IBM
   ILOG CPLEX Optimization Studio 22.1.1.
2. Open `ROernightTrainLogistics/RoOvernightTrainLogistics.sln`, restore the
   NuGet packages, and update the two CPLEX assembly paths in
   `ROernightTrainLogistics/RoOvernightTrainLogistics.csproj` if CPLEX is
   installed in another location.
3. Copy `Input_Nodes.csv`, `Input_Links.csv`, and `Input_Paths.csv` from a
   selected `Results/UploadData/N**D**/data*/` folder to
   `ROernightTrainLogistics/Dataset/`.
4. Copy a desired disruption scenario file to the same directory and name it
   `Input_Disruptions.csv`.
5. Build and run the solution. The program has no command-line arguments and
   writes its output to `TestLog/`.
6. Compare the aggregated output with the relevant worksheet in
   `RobustTrainresults.xlsx`.

## Ongoing Development

This repository is maintained as a replication snapshot for the associated
study. Any post-review corrections will be documented through repository
commits and releases.

## Support

For support in using this software, please open an
[issue](https://github.com/fengtao818/2025.1492.07/issues/new) or contact the
corresponding author listed in [AUTHORS.txt](AUTHORS.txt).
