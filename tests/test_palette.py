"""The two arc colors stay apart for a reader with a color vision deficiency.

Arrowheads are uniform, so a produced and a consumed arc differ by color (and by their direction). grownet's
check, with its simulation (Vienot, Brettel and Mollon 1999): the two colors must stay far apart under
protanopia and deuteranopia, and read as text on white.
"""
import math

import pytest

from foodnet import brand

_RGB_TO_LMS = ((0.31399, 0.63951, 0.04649), (0.15537, 0.75789, 0.08670), (0.01775, 0.10945, 0.87262))
_LMS_TO_RGB = ((5.47221, -4.64196, 0.16963), (-1.12524, 2.29317, -0.16789), (0.02980, -0.19318, 1.16364))
_DICHROMACY = {
    "protanopia": ((0, 1.05118294, -0.05116099), (0, 1, 0), (0, 0, 1)),
    "deuteranopia": ((1, 0, 0), (0.9513092, 0, 0.04866992), (0, 0, 1)),
}
APART = 70.0        # sRGB units between the two colors after simulation (grownet's bar)


def _channels(color: str) -> tuple:
    color = color.lstrip("#")
    return tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))


def _to_linear(value: float) -> float:
    value /= 255
    return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4


def _apply(matrix, vector) -> list:
    return [sum(row[i] * vector[i] for i in range(3)) for row in matrix]


def simulate(color: str, kind: str) -> tuple:
    """`color` as a reader with that dichromacy sees it, as sRGB channels."""
    linear = [_to_linear(c) for c in _channels(color)]
    seen = _apply(_LMS_TO_RGB, _apply(_DICHROMACY[kind], _apply(_RGB_TO_LMS, linear)))
    out = []
    for value in seen:
        value = min(1.0, max(0.0, value))
        value = 12.92 * value if value <= 0.0031308 else 1.055 * value ** (1 / 2.4) - 0.055
        out.append(round(value * 255))
    return tuple(out)


def luminance(color: str) -> float:
    r, g, b = (_to_linear(c) for c in _channels(color))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(color: str, against: str = "#ffffff") -> float:
    high, low = sorted((luminance(color), luminance(against)), reverse=True)
    return (high + 0.05) / (low + 0.05)



@pytest.mark.parametrize("kind", sorted(_DICHROMACY))
def test_produced_and_consumed_stay_apart(kind):
    apart = math.dist(simulate(brand.PRODUCED, kind), simulate(brand.CONSUMED, kind))
    assert apart > APART, f"{brand.PRODUCED} and {brand.CONSUMED} are {apart:.0f} apart under {kind}"


@pytest.mark.parametrize("color", [brand.PRODUCED, brand.CONSUMED])
def test_the_arc_colors_are_readable_as_text_on_white(color):
    assert contrast(color) >= 4.5


def test_foodnet_does_not_reuse_grownets_signal_colors():
    assert {brand.PRODUCED.lower(), brand.CONSUMED.lower()}.isdisjoint({"#1a7f5a", "#c2410c"})


def test_no_genus_color_reads_as_an_arc_color():
    assert brand.GENUS_COLORS and all(brand._away_from_arcs(c) for c in brand.GENUS_COLORS)


@pytest.mark.parametrize("kind", sorted(_DICHROMACY))
@pytest.mark.parametrize("color", ["#ff0000", "#00ff00", "#2160a8", "#b45309"])
def test_the_simulation_lands_on_the_dichromatic_plane(color, kind):
    # a check on the check: without the long or medium cone, red and green come back nearly equal
    red, green, _ = simulate(color, kind)
    assert abs(red - green) <= 12
