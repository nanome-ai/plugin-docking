from nanome.util import Vector3


def get_complex_center(complex):
    """Calculate the center of a complex."""
    inf = float('inf')
    min_pos = Vector3(inf, inf, inf)
    max_pos = Vector3(-inf, -inf, -inf)

    for atom in complex.atoms:
        min_pos.x = min(min_pos.x, atom.position.x)
        min_pos.y = min(min_pos.y, atom.position.y)
        min_pos.z = min(min_pos.z, atom.position.z)
        max_pos.x = max(max_pos.x, atom.position.x)
        max_pos.y = max(max_pos.y, atom.position.y)
        max_pos.z = max(max_pos.z, atom.position.z)

    return (min_pos + max_pos) * 0.5
