# Hybrid Optimization Method for Protein Complex Extraction from PPI Networks Using HPC Platforms

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
- **Infrastructure**: Grid5000, Cytoscape for visualization

## 📊 Results
### Performance Comparison (Yeast)

| Method       | F-measure | Accuracy | Score |
|--------------|-----------|----------|-------|
| OptClust-Hp  | 0.6481    | 0.6552   | 3.2247|
| EWCA         | 0.5622    | 0.5943   | 3.0027|
| COACH        | 0.5022    | 0.5362   | 2.6706|

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
