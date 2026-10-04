import numpy as np
import plotly.graph_objects as go
import periodictable as pt
import apav as ap
import ipywidgets as widgets
import re

from pathlib import Path
from IPython.display import display


# ============================================================
# Elements and maximum atomic charges
# ============================================================

MAIN_ELEMENTS = [
    "Fe", "C", "Mn", "Si",
    "Al", "O", "H", "Ga", "Nb"
]

# Elements shown in the GUI and used for candidate generation.
# N, P, S, Cr, Ni and B are intentionally excluded.
ALL_ELEMENTS = MAIN_ELEMENTS.copy()

DEFAULT_ELEMENTS = MAIN_ELEMENTS.copy()

DEFAULT_MAX_CHARGE = {
    "Fe": 3,
    "C": 2,
    "Mn": 3,
    "Si": 2,
    "Al": 2,
    "O": 2,
    "H": 1,
    "Ga": 2,
    "Nb": 3,
}


# ============================================================
# Isotope functions
# ============================================================

def isotope_records(symbol):

    element = pt.elements.symbol(symbol)
    records = []

    for isotope_number in getattr(element, "isotopes", []):

        isotope = element[isotope_number]
        mass = getattr(isotope, "mass", None)
        abundance = getattr(isotope, "abundance", None)

        if mass is None:
            continue

        if abundance is None:
            abundance = 0.0

        records.append(
            (
                int(isotope_number),
                float(mass),
                float(abundance)
            )
        )

    if not records:
        records.append(
            (
                0,
                float(element.mass),
                100.0
            )
        )

    return records


def most_abundant_isotope(symbol):

    return max(
        isotope_records(symbol),
        key=lambda item: item[2]
    )


def isotope_fingerprint(formula):

    """Return nominal isotope-combination masses for a formula."""

    possible_masses = {0}

    for symbol, atom_count in formula.items():

        isotope_numbers = []

        for isotope_number, isotope_mass, abundance in (
            isotope_records(symbol)
        ):

            if abundance < 0.01:
                continue

            nominal_number = (
                isotope_number
                if isotope_number != 0
                else int(round(isotope_mass))
            )

            isotope_numbers.append(nominal_number)

        for _ in range(atom_count):

            possible_masses = {
                current + isotope_number
                for current in possible_masses
                for isotope_number in isotope_numbers
            }

    return tuple(sorted(possible_masses))


# ============================================================
# Candidate functions
# ============================================================

def formula_mass(formula):

    total = 0.0

    for symbol, count in formula.items():

        _, mass, _ = (
            most_abundant_isotope(symbol)
        )

        total += count * mass

    return total


def add_candidate(
    candidates,
    mass,
    label,
    category,
    abundance=100.0,
    family=None
):

    if mass <= 0:
        return

    candidates.append(
        {
            "mass": float(mass),
            "label": label,
            "category": category,
            "family": family or category,
            "abundance": float(abundance)
        }
    )


def add_molecular_formula(
    candidates,
    formula,
    formula_label,
    category,
    charge=1,
    family=None
):

    base_mass = formula_mass(formula)

    charge_label = (
        "+"
        if charge == 1
        else f"^{charge}+"
    )

    # Main molecular ion
    add_candidate(
        candidates,
        base_mass / charge,
        f"{formula_label}{charge_label}",
        category,
        100.0,
        family=family
    )

    # First-order isotope variants
    for symbol, atom_count in formula.items():

        (
            main_number,
            main_mass,
            main_abundance
        ) = most_abundant_isotope(symbol)

        for (
            isotope_number,
            isotope_mass,
            isotope_abundance
        ) in isotope_records(symbol):

            if isotope_number == main_number:
                continue

            if isotope_abundance < 0.01:
                continue

            variant_mass = (
                base_mass
                + isotope_mass
                - main_mass
            )

            relative_abundance = (
                100.0
                * atom_count
                * isotope_abundance
                / main_abundance
            )

            add_candidate(
                candidates,
                variant_mass / charge,
                (
                    f"{formula_label}{charge_label} "
                    f"({isotope_number}{symbol})"
                ),
                category,
                relative_abundance,
                family=family
            )


def parse_manual_candidates(text):

    result = []

    for item in text.split(";"):

        item = item.strip()

        if not item:
            continue

        if ":" not in item:
            raise ValueError(
                "Manual peaks must use label:mass"
            )

        label, mass_text = item.split(":", 1)

        result.append(
            {
                "mass": float(mass_text.strip()),
                "label": label.strip(),
                "category": "Manual candidates",
                "abundance": 100.0
            }
        )

    return result


def extract_charge(label):

    core_label = label.split(" ", 1)[0]

    match = re.search(r"\^(\d+)\+$", core_label)

    if match:
        return int(match.group(1))

    if core_label.endswith("+"):
        return 1

    return ""


# ============================================================
# Individual molecular-ion options
# ============================================================

def build_molecular_options():

    options = []

    def add(
        key,
        label,
        formula,
        category,
        charge=1,
        default=False
    ):

        fingerprint_text = ", ".join(
            str(number)
            for number in isotope_fingerprint(formula)
        )

        options.append(
            {
                "key": key,
                "label": label,
                "gui_label": f"{label} ({fingerprint_text})",
                "formula": formula,
                "category": category,
                "charge": charge,
                "default": default
            }
        )

    molecular_elements = [
        symbol
        for symbol in ALL_ELEMENTS
        if symbol not in ["C", "H", "O"]
    ]

    # H2+
    add(
        "H2",
        "H2+",
        {"H": 2},
        "Hydrogen molecular ions",
        default=True
    )

    # Carbon clusters
    for count in [2, 3, 4, 5, 6]:

        add(
            f"C{count}",
            f"C{count}+",
            {"C": count},
            "Carbon clusters",
            default=count <= 3
        )

    # Doubly charged carbon clusters
    for count in [2, 3]:

        add(
            f"C{count}_2plus",
            f"C{count}^2+",
            {"C": count},
            "Doubly charged carbon clusters",
            charge=2
        )

    # Oxygen clusters
    for count in [2, 3]:

        add(
            f"O{count}",
            f"O{count}+",
            {"O": count},
            "Oxygen molecular ions",
            default=count == 2
        )

    for symbol in molecular_elements:

        # Hydride
        add(
            f"{symbol}H",
            f"{symbol}H+",
            {
                symbol: 1,
                "H": 1
            },
            "Metal hydrides"
        )

        # Carbides
        for carbon_count in [1, 2]:

            add(
                f"{symbol}C{carbon_count}",
                f"{symbol}C{carbon_count}+",
                {
                    symbol: 1,
                    "C": carbon_count
                },
                "Metal-carbides",
                default=(
                    symbol == "Fe"
                    and carbon_count == 1
                )
            )

        # Oxides and doubly charged oxides
        for oxygen_count in [1, 2]:

            add(
                f"{symbol}O{oxygen_count}",
                f"{symbol}O{oxygen_count}+",
                {
                    symbol: 1,
                    "O": oxygen_count
                },
                "Metal-oxides",
                default=(
                    symbol == "Fe"
                    and oxygen_count == 1
                )
            )

            add(
                f"{symbol}O{oxygen_count}_2plus",
                f"{symbol}O{oxygen_count}^2+",
                {
                    symbol: 1,
                    "O": oxygen_count
                },
                "Doubly charged metal oxides",
                charge=2
            )

        # Hydroxide
        add(
            f"{symbol}OH",
            f"{symbol}OH+",
            {
                symbol: 1,
                "O": 1,
                "H": 1
            },
            "Metal hydroxides"
        )

    # Carbon-hydrogen ions
    for label, formula in [
        ("CH", {"C": 1, "H": 1}),
        ("CH2", {"C": 1, "H": 2}),
        ("C2H", {"C": 2, "H": 1})
    ]:

        add(
            label,
            f"{label}+",
            formula,
            "Carbon-hydrogen ions"
        )

    # Mixed-metal ions
    for index, first in enumerate(molecular_elements):

        for second in molecular_elements[index + 1:]:

            label = f"{first}{second}"

            add(
                label,
                f"{label}+",
                {
                    first: 1,
                    second: 1
                },
                "Mixed-metal ions"
            )

    return options


MOLECULAR_OPTIONS = build_molecular_options()


# ============================================================
# Generate candidates
# ============================================================

def generate_candidates(
    selected_elements,
    charge_limits,
    include_atomic,
    selected_molecular_options,
    max_cluster_size
):

    candidates = []

    # Atomic ions
    if include_atomic:

        for symbol in selected_elements:

            maximum_charge = charge_limits[symbol]

            for charge in range(
                1,
                maximum_charge + 1
            ):

                for (
                    isotope_number,
                    exact_mass,
                    abundance
                ) in isotope_records(symbol):

                    if abundance <= 0:
                        continue

                    isotope_label = (
                        symbol
                        if isotope_number == 0
                        else f"{isotope_number}{symbol}"
                    )

                    charge_label = (
                        "+"
                        if charge == 1
                        else f"^{charge}+"
                    )

                    add_candidate(
                        candidates,
                        exact_mass / charge,
                        f"{isotope_label}{charge_label}",
                        f"Atomic ions, charge {charge}+",
                        abundance,
                        family=f"{symbol}{charge_label}"
                    )

    # Individual molecular ions
    for option in selected_molecular_options:

        formula = option["formula"]

        # Skip if an element is not selected
        if not all(
            symbol in selected_elements
            for symbol in formula
        ):
            continue

        # Limit carbon cluster size
        if formula.get("C", 0) > max_cluster_size:
            continue

        label = option["label"]

        if label.endswith("^2+"):
            formula_label = label[:-3]
        elif label.endswith("+"):
            formula_label = label[:-1]
        else:
            formula_label = label

        add_molecular_formula(
            candidates,
            formula,
            formula_label,
            option["category"],
            option["charge"],
            family=option["label"]
        )

    # Remove duplicates
    unique = {}

    for candidate in candidates:

        key = (
            round(candidate["mass"], 5),
            candidate["label"]
        )

        unique[key] = candidate

    candidates = list(unique.values())

    candidates.sort(
        key=lambda item: item["mass"]
    )

    return candidates


# ============================================================
# Plot
# ============================================================

def make_spectrum_plot(
    roi,
    candidates,
    lower_mass,
    upper_mass,
    bin_width,
    logarithmic_y,
    show_labels
):

    mass, counts = roi.mass_histogram(
        bin_width=bin_width,
        lower=lower_mass,
        upper=upper_mass,
        multiplicity="all"
    )

    positive = counts > 0

    if not np.any(positive):
        raise ValueError(
            "No positive counts found in this mass range."
        )

    mass = mass[positive]
    counts = counts[positive]

    y_min = max(
        1,
        float(np.min(counts))
    )

    y_max = float(
        np.max(counts)
    )

    fig = go.Figure()

    # Measured spectrum
    fig.add_trace(
        go.Scattergl(
            x=mass,
            y=counts,
            mode="lines",
            line=dict(
                color="black",
                width=1
            ),
            name="Measured spectrum",
            hovertemplate=(
                "Mass: %{x:.5f} Da<br>"
                "Counts: %{y:.0f}"
                "<extra></extra>"
            )
        )
    )

    dashes = [
        "solid",
        "dot",
        "dash",
        "dashdot",
        "longdash",
        "longdashdot"
    ]

    families = sorted(
        set(
            item.get("family", item["category"])
            for item in candidates
        )
    )

    styles = {}

    for index, family in enumerate(families):

        hue = (index * 360.0 / max(len(families), 1)) % 360.0

        styles[family] = {
            "color": f"hsl({hue:.1f}, 75%, 45%)",
            "dash": dashes[
                (index // 18)
                % len(dashes)
            ]
        }

    shown_in_legend = set()

    # Candidate guide sticks
    for candidate in candidates:

        candidate_mass = candidate["mass"]

        if not (
            lower_mass
            <= candidate_mass
            <= upper_mass
        ):
            continue

        abundance = min(
            max(
                candidate.get(
                    "abundance",
                    100.0
                ),
                0.0
            ),
            100.0
        )

        line_top = (
            y_min
            + abundance / 100.0
            * (y_max - y_min)
        )

        family = candidate.get(
            "family",
            candidate["category"]
        )
        style = styles[family]

        description = (
            f"<b>{candidate['label']}</b><br>"
            f"Family: {family}<br>"
            f"Category: {candidate['category']}<br>"
            f"Mass-to-charge: "
            f"{candidate_mass:.5f} Da<br>"
            f"Natural abundance: "
            f"{abundance:.3f}%"
        )

        show_legend = (
            family not in shown_in_legend
        )

        shown_in_legend.add(family)

        fig.add_trace(
            go.Scatter(
                x=[candidate_mass, candidate_mass],
                y=[y_min, line_top],
                mode="lines",
                line=dict(
                    color=style["color"],
                    width=1.8,
                    dash=style["dash"]
                ),
                opacity=0.95,
                name=family,
                legendgroup=family,
                showlegend=show_legend,
                hovertemplate=(
                    description
                    + "<extra></extra>"
                )
            )
        )

    # Labels
    if show_labels:

        label_x = []
        label_y = []
        label_text = []

        for candidate in candidates:

            candidate_mass = candidate["mass"]

            if (
                lower_mass
                <= candidate_mass
                <= upper_mass
            ):

                label_x.append(candidate_mass)
                label_y.append(y_max)
                label_text.append(candidate["label"])

        fig.add_trace(
            go.Scattergl(
                x=label_x,
                y=label_y,
                mode="text",
                text=label_text,
                textposition="top center",
                textfont=dict(size=9),
                name="Candidate labels",
                showlegend=False,
                hoverinfo="skip"
            )
        )

    fig.update_layout(
        title=(
            "APT Candidate List Maker - "
            "Mass Spectrum"
        ),
        template="plotly_white",
        width=1800,
        height=1050,
        dragmode="pan",
        hovermode="closest",
        margin=dict(
            t=190,
            r=330,
            b=80,
            l=90
        ),
        legend=dict(
            title="Ion family",
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top"
        )
    )

    fig.update_xaxes(
        title="Mass-to-charge ratio (Da)",
        range=[lower_mass, upper_mass],
        rangeslider_visible=True
    )

    fig.update_yaxes(
        title="Counts",
        type=(
            "log"
            if logarithmic_y
            else "linear"
        )
    )

    return fig
