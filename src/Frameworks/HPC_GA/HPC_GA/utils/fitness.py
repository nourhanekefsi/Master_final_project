def normalize(fitness_values: list) -> list:
    max_val = max(fitness_values)
    return [f / max_val for f in fitness_values]