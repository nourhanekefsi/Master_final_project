import logging
import multiprocessing
import sys
import os
from typing import List, Set
import ray

from HPC_GA.parallel import IslandModel
from data_loader import (
    load_ppi, load_complexes,
    random_initialization, load_initialization_complexes
)
from genetic_algorithm import AlgorithmGA, ComplexCrossover, TabuMutator_MasterSlave
from population import MyPopulation


# Logger global
logger = logging.getLogger("MainLogger")

def setup_logger(log_file: str, err_file: str):
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    os.makedirs(os.path.dirname(err_file), exist_ok=True)

    logger.setLevel(logging.INFO)
    logger.propagate = False  # Empêche la duplication via le root logger
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Clear existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    fh = logging.FileHandler(log_file, mode='a')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    eh = logging.FileHandler(err_file, mode='a')
    eh.setFormatter(formatter)
    logger.addHandler(eh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)


def save_results(solution: List[Set[str]], output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for cluster in solution:
            f.write(" ".join(cluster) + "\n")

def main():
    try:
        # Setup logging
        setup_logger("kefsi_mekhazni_workspace/results/STRING_humain/string_humain.log",
                     "kefsi_mekhazni_workspace/results/STRING_humain/string_humain.err")

        if ray.is_initialized():
            ray.shutdown()

        logger.info("Initialisation de Ray vers le cluster")
        ray.init(address="auto", logging_level=logging.INFO, ignore_reinit_error=True)

        resources = ray.cluster_resources()
        if resources.get("CPU", 0) <= 1:
            raise RuntimeError("Pas assez de ressources CPU disponibles.")
        logger.info(f"Cluster prêt. Ressources : {resources}")

        BASE_DIR = os.path.expanduser("~/kefsi_mekhazni_workspace")
        for sub in ["data", "results"]:
            os.makedirs(os.path.join(BASE_DIR, sub), exist_ok=True)

        PPI_PATH = os.path.join(BASE_DIR, "data/ppi_networks/weighted_STRING_humain.txt")
        REF_PATH = os.path.join(BASE_DIR, "data/golden_complexes/STRING_humain.txt")
        INIT_PATH = os.path.join(BASE_DIR, "data/initialization_complexes/detected_complexes_STRING_humain.txt")
        RESULT_PATH = os.path.join(BASE_DIR, "results//STRING_humain/best_solution.txt")

        for p in [PPI_PATH, REF_PATH, INIT_PATH]:
            if not os.path.exists(p):
                raise FileNotFoundError(f"Fichier introuvable : {p}")

        ppi = load_ppi(PPI_PATH)
        logger.info(f"Le réseau PPI a été chargé avec succès. Taille: {len(ppi)}")

        references = load_complexes(REF_PATH)
        logger.info(f"Les complexes de référence ont été chargés. Nombre: {len(references)}")

        population1 = load_initialization_complexes(INIT_PATH, ppi)
        logger.info(f"Population initiale chargée. Taille: {len(population1)}")

        logger.info("Debut de la génération de la population aléatoire")
        population2 = random_initialization(references, ppi, population_size=26,
                                            max_complexes=len(references)*10, max_size=25)
        logger.info(f"Population aléatoire générée. Taille: {len(population2)}")

        total_population = population1 + population2
        population = MyPopulation(total_population)
        logger.info(f"Population totale créée. Taille: {len(population.individuals)}")

        nb_islands = 3
        crossover = ComplexCrossover(ppi)
        # Calcul du nombre de CPU total dans le cluster
        total_cluster_cpus = sum(node["Resources"].get("CPU", 0) for node in ray.nodes() if node["Alive"])
        logger.info(f"Nombre total de CPU dans le cluster : {total_cluster_cpus}")
        
        if total_cluster_cpus <= 6:
            raise RuntimeError("Pas assez de CPU disponibles dans le cluster après réservation.")

        nb_processors = int(max(1, (total_cluster_cpus - 6) // nb_islands))
        logger.info(f"Nombre de processeurs attribués à chaque mutator : {nb_processors}")

        mutator = TabuMutator_MasterSlave(
            graph=ppi,
            tabu_tenure=20,
            tabu_iterations=100,
            rate=1,
            threshold_ratio=600,
            nb_processors=nb_processors
     )


        ga_config = {
            "references": references,
            "ppi_name": "STRING_humain",
            "graph": ppi,
            "population": population,
            "crossover": crossover,
            "mutator": mutator,
            "selection": population.tournament_selection,
            "update_population": population._update_population,
            "max_generations": 50,
            "nb_parents": 2,
            "update_type": "replace"
        }

        parallel_config = {
            "n_islands": nb_islands,
            "migration_interval": 5,
            "migration_size": 3,
            "migration_topology": "ring",
            "migration_config": {
                "strategy": "worst"
            }
        }

        ga = IslandModel(AlgorithmGA, ga_config, parallel_config)
        best_solution = ga.run(generations=ga_config["max_generations"])
        save_results(best_solution.genes, RESULT_PATH)

        logger.info(f"Fitness de la meilleure solution trouvée: {best_solution.evaluate()}")
        logger.info(f"Résultats sauvegardés dans : {RESULT_PATH}")

    except Exception as e:
        logger.exception("Erreur détectée lors de l'exécution principale")
        raise
    finally:
        if ray.is_initialized():
            ray.shutdown()
        logger.info("Ray arrêté. Terminé.")

if __name__ == "__main__":
    main()
