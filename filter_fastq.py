

def filter_gc_content(seqs: Dict[str, Tuple[str, str]],
                      gc_bounds: Union[float, Tuple[Union[float, int], Union[float, int]]] = (0, 100)) -> Dict[str, float]:
    """
    Filters sequences based on GC-content.

    Parameters:
    -----------
    seqs (Dict[str, Tuple[str, str]]):
        A dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    gc_bounds (Union[float, Tuple[Union[float, int], Union[float, int]]]):
        GC-content filtering bounds.
        - Default is 100 (no filtering).
        - If a single float is provided, it's considered as the upper bound.
        - If an integer is provided, it's considered as the upper bound.
        - If a tuple of two numbers is provided, they are treated as the lower and upper bounds.

    Returns:
    --------
    Dict[str, float]:
        A filtered dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    """
    filtered_seqs = {}

    if isinstance(gc_bounds, (float, int)):
        lower_bound = 0
        upper_bound = gc_bounds
    elif isinstance(gc_bounds, tuple) and len(gc_bounds) == 2:
        lower_bound, upper_bound = gc_bounds
    else:
        raise ValueError("gc_bounds should be a single float/int or a tuple of two floats/ints.")

    for key, value in seqs.items():
        gc_count = 0
        not_gc_count = 0

        for base in value[0]:
            if base in ('G', 'C'):
                gc_count += 1
            else:
                not_gc_count += 1

        gc_content = gc_count / (gc_count + not_gc_count) * 100

        if lower_bound <= gc_content <= upper_bound:
            filtered_seqs[key] = value

    return filtered_seqs

def filter_length(seqs: Dict[str, Tuple[str, str]],
                  length_bounds: Union[int, Tuple[int, int]] = (0, 2**32)) -> Dict[str, int]:
    """
    Filters sequences based on sequence length.

    Parameters:
    -----------
    seqs (Dict[str, Tuple[str, str]]):
        A dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    length_bounds (Union[int, Tuple[int, int]]):
        Length filtering bounds.
        - Default is (0, 2^32) (no filtering).
        - If a single integer is provided, it's considered as the upper bound.
        - If a tuple of two integers is provided, they are treated as the lower and upper bounds.

    Returns:
    --------
    Dict[str, int]:
        A filtered dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.

    """
    filtered_seqs = {}

    if isinstance(length_bounds, int):
        lower_bound = 0
        upper_bound = length_bounds
    elif isinstance(length_bounds, tuple) and len(length_bounds) == 2:
        lower_bound, upper_bound = length_bounds
    else:
        raise ValueError("length_bounds should be a single integer or a tuple of two integers.")

    for key, value in seqs.items():
        seq_length = len(value[0])

        if lower_bound <= seq_length <= upper_bound:
            filtered_seqs[key] = value

    return filtered_seqs
