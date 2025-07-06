# Hybrid Optimization Method for Protein Complex Extraction from PPI Networks Using HPC Platforms

## 📚 Table of Contents
- [Overview](#-overview)
- [Objectives](#-objectives)
- [Methodology](#-methodology)
- [Technologies Used](#-technologies-used)
- [Results](#-results)
- [Repository Structure](#-repository-structure)

## 📌 Overview
This project detects protein complexes from large-scale protein-protein interaction (PPI) networks using a hybrid optimization approach combining:
- PPI network weighting with biological data
- Multi-objective fitness evaluation
- Memetic algorithm (Genetic Algorithm + Tabu Search)
- Parallel computing via island model and master-slave architecture

## 🎯 Objectives
- Automate complex detection in large PPI networks
- Address limitations of existing methods (noisy PPIs, complex structure variability)
- Improve accuracy through hybrid metaheuristics
- Reduce computation time via parallel execution
- Validation using real PPI datasets and gold-standard complexes

## 🧪 Methodology
### Hybrid Approach
- **Genetic Algorithm** for exploration
- **Tabu Search** for local optimization

### Data Processing
1. **Data Sources**:
   - PPI databases: BioGRID (yeast/human), STRING (yeast/human), DIP (yeast)
   - Gold standards: CORUM, Hu.Map, SGD, Complex Portal
   - Biological data: GO slim, subcellular localization, gene expression

2. **Encoding**:
   - Complex = Protein set (e.g., C₁ = {P₁,P₂,P₃})
   - Solution = Set of complexes

3. **Fitness Function**:
F(C) = Cohesion + Density + Isolation - Boundary + Modularity

### Parallel Implementation
- **Island model** with Ray/ThreadPool
- Each island contains:
   - Tabu workers (complex optimization)
   - Master-slave structure for metric computation

## 🛠️ Technologies
- **Python** with NumPy, Pandas, Matplotlib
- **HPC Frameworks**:
- `hpc-ga`: Supports islands, cellular, master-slave models
- `hpc-tabu`: Implements Crainic et al. taxonomy
- ray, multiprocessing
- **Infrastructure**: Grid5000, Cytoscape for visualization

## 📊 Results
### Metrics – BIOGRID (Human)
| Method       | F-Measure | Accuracy | Jaccard | CR    | MMR   | Total Score |
|--------------|-----------|----------|---------|-------|-------|--------------|
| MCL          | 0.4903    | 0.4965   | 0.2761  | 0.8136| 0.3592| 2.4357       |
| MCODE        | 0.3883    | 0.3954   | 0.2522  | 0.7598| 0.3642| 2.1599       |
| ClusterONE   | 0.2989    | 0.3191   | 0.1425  | 0.9738| 0.2362| 1.9704       |
| COACH        | 0.5022    | 0.5362   | 0.2544  | 0.9770| 0.4008| 2.6706       |
| SE-DMTG      | 0.4404    | 0.4788   | 0.3932  | 0.5466| 0.5548| 2.4138       |
| EWCA         | 0.5622    | 0.5943   | 0.3810  | 0.9396| 0.5256| 3.0027       |
| MPC-C        | 0.4910    | 0.5310   | 0.3104  | 0.9482| 0.4723| 2.7529       |
| OptClust-H   | 0.6300    | 0.6400   | 0.4000  | 0.5400| 0.9300| 3.1600       |
| OptClust-Hp  | 0.6481    | 0.6552   | 0.4153  | 0.5610| 0.9451| 3.2247       |

### Metrics – STRING (Human)
| Method       | F-Measure | Accuracy | Jaccard | CR    | MMR   | Total Score |
|--------------|-----------|----------|---------|-------|-------|--------------|
| MCL          | 0.3656    | 0.4235   | 0.1735  | 0.8565| 0.2357| 2.0548       |
| MCODE        | 0.4026    | 0.4417   | 0.1982  | 0.8869| 0.2893| 2.2187       |
| ClusterONE   | 0.2524    | 0.2922   | 0.1170  | 0.9587| 0.1998| 1.8201       |
| COACH        | 0.4063    | 0.4792   | 0.1851  | 0.9838| 0.3086| 2.3630       |
| SE-DMTG      | 0.3742    | 0.4351   | 0.2112  | 0.6906| 0.3232| 2.0343       |
| EWCA         | 0.4042    | 0.4876   | 0.2167  | 0.9878| 0.3551| 2.4514       |
| MPC-C        | 0.3640    | 0.4450   | 0.2040  | 0.9513| 0.3234| 2.2877       |
| OptClust-H   | 0.4300    | 0.5000   | 0.2200  | 0.3300| 0.9900| 2.4700       |
| OptClust-Hp  | 0.4440    | 0.5138   | 0.2351  | 0.3557| 0.9931| 2.5417       |

### Metrics – STRING (Yeast) 

| Méthode      | F-measure | Accuracy   | Jaccard | CR     | MMR    | Total score |
|--------------|----------|--------|---------|--------|--------|--------------|
| MCL          | 0.3398   | 0.4340 | 0.1655  | 0.9205 | 0.2168 | 2.0767       |
| MCODE        | 0.3640   | 0.4381 | 0.1808  | 0.9702 | 0.2483 | 2.2014       |
| ClusterONE   | 0.2234   | 0.2922 | 0.0985  | 0.9139 | 0.1797 | 1.7076       |
| COACH        | 0.3234   | 0.4304 | 0.1557  | 0.9603 | 0.2629 | 2.1327       |
| SE-DMTG      | 0.3052   | 0.3784 | 0.1550  | 0.6490 | 0.2227 | 1.7103       |
| EWCA         | 0.2991   | 0.4161 | 0.1830  | 1.0000 | 0.3004 | 2.1986       |
| MPC-C        | 0.3575   | 0.4527 | 0.1952  | 0.9404 | 0.2833 | 2.2292       |
| OptClust-H   | 0.32     | 0.42   | 0.19    | 0.28   | 0.99   | 2.21         |
| OptClust-Hp  | 0.3428   | 0.4562 | 0.2067  | 0.2974 | 0.9961 | 2.2992       |


### Metrics – BIOGRID (Yeast)
| Method       | F-Measure | Accuracy | Jaccard | CR    | MMR   | Total Score |
|--------------|-----------|----------|---------|-------|-------|--------------|
| MCL          | 0.3838    | 0.4293   | 0.1888  | 0.8436| 0.2414| 2.0869       |
| MCODE        | 0.3294    | 0.3523   | 0.1607  | 0.8601| 0.2417| 1.9442       |
| ClusterONE   | 0.2369    | 0.2841   | 0.1038  | 0.9753| 0.1742| 1.7743       |
| COACH        | 0.4087    | 0.4847   | 0.1954  | 0.9877| 0.3151| 2.3915       |
| SE-DMTG      | 0.3177    | 0.4044   | 0.1641  | 0.7984| 0.2717| 1.9562       |
| EWCA         | 0.4314    | 0.5067   | 0.2341  | 0.9835| 0.3446| 2.5004       |
| MPC-C        | 0.4179    | 0.4879   | 0.2376  | 0.9712| 0.3556| 2.4701       |
| OptClust-H   | 0.4500    | 0.5100   | 0.2300  | 0.3100| 0.9900| 2.5100       |
| OptClust-Hp  | 0.4682    | 0.5283   | 0.2476  | 0.3260| 0.9934| 2.5635       |

### Metrics – DIP (Yeast)
| Method       | F-Measure | Accuracy | Jaccard | CR    | MMR   | Total Score |
|--------------|-----------|----------|---------|-------|-------|--------------|
| MCL          | 0.3400    | 0.3800   | 0.1600  | 0.8200| 0.2000| 1.9300       |
| MCODE        | 0.3300    | 0.3500   | 0.1800  | 0.8100| 0.2500| 1.9400       |
| ClusterONE   | 0.1900    | 0.2000   | 0.0700  | 1.0000| 0.1200| 1.6500       |
| COACH        | 0.3500    | 0.4400   | 0.1600  | 0.9800| 0.2400| 2.1900       |
| SE-DMTG      | 0.1700    | 0.2800   | 0.1700  | 0.2800| 0.7400| 1.6500       |
| EWCA         | 0.3700    | 0.4500   | 0.2300  | 0.3200| 0.9800| 2.3500       |
| MPC-C        | 0.3400    | 0.4200   | 0.2200  | 0.3200| 0.9500| 2.2700       |
| OptClust-H   | 0.2300    | 0.3300   | 0.1900  | 0.2900| 0.9400| 2.0000       |
| OptClust-Hp  | 0.2507    | 0.3592   | 0.2041  | 0.3115| 0.9552| 2.0807       |


### Execution Time
![Parallel vs Sequential Time Comparison](https://github.com/nourhanekefsi/Master_final_project/blob/main/execution_time.png)

## 📂 Repository Structure
/data

      /clean

            /complexes

            /interactions

            /weighted-networks

/notebooks

     /sequential
/src

     /frameworks

        HPC_ga/

        HPC_tabu/

     /parallel

      /sequential

requirements.txt
