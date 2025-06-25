from math import sqrt
from typing import Dict, List, Set
import numpy as np
import ray
import logging
logger = logging.getLogger(__name__)

def overlap_score(pred: List[str], ref: List[str]) -> float:
    return len(set(pred) & set(ref)) / sqrt(len(pred) * len(ref)) if pred and ref else 0.0

def jaccard_index(pred: List[str], ref: List[str]) -> float:
    union = set(pred) | set(ref)
    return len(set(pred) & set(ref)) / len(union) if union else 0.0


@ray.remote
def compute_ppv(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    numerator = sum(max(len(set(pred) & set(ref)) for ref in real_complexes) for pred in predicted_complexes)
    denominator = sum(len(pred) for pred in predicted_complexes)
    return numerator / denominator if denominator > 0 else 0.0

@ray.remote
def compute_sn(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    numerator = sum(max(len(set(pred) & set(ref)) for pred in predicted_complexes) for ref in real_complexes)
    denominator = sum(len(ref) for ref in real_complexes)
    return numerator / denominator if denominator > 0 else 0.0

@ray.remote
def compute_f_measure(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    ppv = ray.get(compute_ppv.remote(predicted_complexes, real_complexes))
    sn = ray.get(compute_sn.remote(predicted_complexes, real_complexes))
    return (2 * ppv * sn) / (ppv + sn) if (ppv + sn) > 0 else 0.0

@ray.remote
def compute_accuracy(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    ppv = ray.get(compute_ppv.remote(predicted_complexes, real_complexes))
    sn = ray.get(compute_sn.remote(predicted_complexes, real_complexes))
    return sqrt(sn * ppv) if ppv > 0 and sn > 0 else 0.0

@ray.remote
def compute_mmr(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    if not predicted_complexes:
        return 0.0
    return np.mean([max(overlap_score(pred, ref) for ref in real_complexes) for pred in predicted_complexes])

@ray.remote
def compute_jaccard(predicted_complexes: List[List[str]], real_complexes: List[List[str]]) -> float:
    if not predicted_complexes:
        return 0.0
    return np.mean([max(jaccard_index(pred, ref) for ref in real_complexes) for pred in predicted_complexes])

@ray.remote
def compute_covered_rate(predicted_complexes: List[List[str]], real_complexes: List[List[str]], threshold: float = 0.2) -> float:
    if not real_complexes:
        return 0.0
    matched_reals = sum(
        any(overlap_score(pred, ref) >= threshold for pred in predicted_complexes)
        for ref in real_complexes
    )
    return matched_reals / len(real_complexes)


class MasterSlaveMetrics:
    def __init__(self, predicted_complexes: List[List[str]], references: List[List[str]], n_workers: int = None):
        self.workers = [
            compute_accuracy,
            compute_f_measure,
            compute_jaccard,
            compute_mmr,
            compute_ppv,
            compute_sn,
            compute_covered_rate
        ]
        self.metric_names = [
            "Accuracy",
            "F-mesure",
            "Jaccard",
            "MMR",
            "PPV",
            "Recall (Sn)",
            "Covered Rate"
        ]
        self.references = references
        self.predicted_complexes = predicted_complexes

    def compute_all_metrics(self) -> Dict[str, float]:
        if not self.predicted_complexes or not self.references:
            return {name: 0.0 for name in self.metric_names}
    
        results = {}
        for metric_function, metric_name in zip(self.workers, self.metric_names):
            try:
                results[metric_name] = ray.get(metric_function.remote(self.predicted_complexes, self.references))
            except Exception as e:
                logger.error(f"Erreur calcul métrique {metric_name}: {str(e)}")
                 
                results[metric_name] = 0.0
            
        results["Score Total"] = sum([results["Accuracy"], results["F-mesure"], results["Jaccard"], results["MMR"], results["Covered Rate"]]) 
        return results
