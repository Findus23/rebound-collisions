def clamp(n: float, smallest: float, largest: float) -> float:
    assert smallest < largest
    return max(smallest, min(n, largest))

