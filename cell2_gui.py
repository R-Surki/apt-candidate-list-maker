# ============================================================
# GUI widgets
# ============================================================

epos_widget = widgets.Text(
    value=(
        "/Users/rsurkial23/Desktop/Experimental data/"
        "Rooh/APT/KTH/"
        "865-CR680C-1s- carbide- low allyed austenite "
        "/Alex/865.EPOS"
    ),
    description="EPOS file:",
    layout=widgets.Layout(width="100%")
)

lower_widget = widgets.FloatText(
    value=0.0,
    description="Lower m/z:"
)

upper_widget = widgets.FloatText(
    value=200.0,
    description="Upper m/z:"
)

bin_width_widget = widgets.FloatText(
    value=0.01,
    description="Bin width:"
)

cluster_size_widget = widgets.IntSlider(
    value=4,
    min=2,
    max=8,
    step=1,
    description="Max C cluster:"
)

atomic_widget = widgets.Checkbox(
    value=True,
    description="Atomic ions"
)

log_widget = widgets.Checkbox(
    value=True,
    description="Logarithmic y-axis"
)

labels_widget = widgets.Checkbox(
    value=True,
    description="Show labels"
)

manual_widget = widgets.Textarea(
    value="",
    description="Manual peaks:",
    placeholder=(
        "Example: unknown-1: 37.123; "
        "possible-ion: 44.987"
    ),
    layout=widgets.Layout(
        width="100%",
        height="80px"
    )
)

save_widget = widgets.Text(
    value=(
        "/Users/rsurkial23/Desktop/"
        "individual_candidate_spectrum.html"
    ),
    description="Save HTML:",
    layout=widgets.Layout(width="100%")
)

candidate_save_widget = widgets.Text(
    value=(
        "/Users/rsurkial23/Desktop/"
        "candidate_list.txt"
    ),
    description="Save candidates:",
    layout=widgets.Layout(width="100%")
)

button = widgets.Button(
    description="Generate spectrum",
    button_style="success",
    icon="line-chart"
)

output = widgets.Output()

# Keep the EPOS data in notebook memory so changing candidate selections
# does not reload the file every time the button is pressed.
cached_roi = None
cached_epos_signature = None


# Element and charge widgets
element_widgets = {}
charge_widgets = {}

for symbol in ALL_ELEMENTS:

    element_widgets[symbol] = widgets.Checkbox(
        value=(
            symbol in DEFAULT_ELEMENTS
        ),
        description=symbol,
        indent=False,
        layout=widgets.Layout(width="65px")
    )

    charge_widgets[symbol] = widgets.BoundedIntText(
        value=DEFAULT_MAX_CHARGE[symbol],
        min=1,
        max=6,
        step=1,
        description="max:",
        layout=widgets.Layout(width="125px")
    )


# Molecular-ion widgets
molecular_widgets = {}

for option in MOLECULAR_OPTIONS:

    molecular_widgets[option["key"]] = (
        widgets.Checkbox(
            value=option["default"],
            description=option["gui_label"],
            indent=False,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="320px")
        )
    )


# Group molecular-ion widgets
molecular_categories = {}

for option in MOLECULAR_OPTIONS:

    category = option["category"]

    molecular_categories.setdefault(
        category,
        []
    ).append(
        molecular_widgets[option["key"]]
    )


molecular_boxes = []

for category, category_widgets in (
    molecular_categories.items()
):

    rows = []

    for start in range(
        0,
        len(category_widgets),
        3
    ):

        rows.append(
            widgets.HBox(
                category_widgets[
                    start:start + 5
                ]
            )
        )

    molecular_boxes.append(
        widgets.VBox(
            [
                widgets.HTML(
                    f"<b>{category}</b>"
                ),
                *rows
            ]
        )
    )


# ============================================================
# Button callback
# ============================================================

def generate_button_clicked(_):

    global cached_roi
    global cached_epos_signature

    with output:

        output.clear_output()

        try:

            epos_path = Path(
                epos_widget.value.strip()
            )

            if not epos_path.is_file():

                raise FileNotFoundError(
                    f"EPOS file not found:\n{epos_path}"
                )

            selected_elements = [
                symbol
                for symbol in ALL_ELEMENTS
                if element_widgets[symbol].value
            ]

            if not selected_elements:

                raise ValueError(
                    "Select at least one element."
                )

            charge_limits = {
                symbol: charge_widgets[symbol].value
                for symbol in selected_elements
            }

            selected_options = [
                option
                for option in MOLECULAR_OPTIONS
                if molecular_widgets[
                    option["key"]
                ].value
            ]

            if (
                lower_widget.value
                >= upper_widget.value
            ):

                raise ValueError(
                    "Lower mass must be smaller "
                    "than upper mass."
                )

            if bin_width_widget.value <= 0:

                raise ValueError(
                    "Bin width must be greater than zero."
                )

            # Reload only when the path or the file itself has changed.
            file_stat = epos_path.stat()
            epos_signature = (
                str(epos_path.resolve()),
                file_stat.st_mtime_ns,
                file_stat.st_size
            )

            if (
                cached_roi is None
                or cached_epos_signature != epos_signature
            ):

                print("Loading EPOS file...")

                cached_roi = ap.Roi.from_epos(
                    str(epos_path)
                )

                cached_epos_signature = epos_signature

            else:

                print("Reusing EPOS data already loaded in memory...")

            roi = cached_roi

            candidates = generate_candidates(
                selected_elements=selected_elements,
                charge_limits=charge_limits,
                include_atomic=atomic_widget.value,
                selected_molecular_options=(
                    selected_options
                ),
                max_cluster_size=(
                    cluster_size_widget.value
                )
            )

            candidates.extend(
                parse_manual_candidates(
                    manual_widget.value
                )
            )

            candidates.sort(
                key=lambda item: item["mass"]
            )

            candidate_save_path = (
                candidate_save_widget.value.strip()
            )

            if candidate_save_path:

                candidate_lines = [
                    "Ion name\tCharge\tm/z (Da)\t"
                    "Family\tCategory\tAbundance (%)"
                ]

                for candidate in candidates:

                    candidate_lines.append(
                        "\t".join(
                            [
                                candidate["label"],
                                str(
                                    extract_charge(
                                        candidate["label"]
                                    )
                                ),
                                f"{candidate['mass']:.8f}",
                                candidate.get(
                                    "family",
                                    candidate["category"]
                                ),
                                candidate["category"],
                                f"{candidate.get('abundance', 100.0):.6f}"
                            ]
                        )
                    )

                Path(candidate_save_path).write_text(
                    "\n".join(candidate_lines),
                    encoding="utf-8"
                )

                print(
                    f"Saved candidate list to:\n"
                    f"{candidate_save_path}"
                )

            print(
                f"Loaded {roi.counts:,} ions"
            )

            print(
                f"Selected elements: "
                f"{selected_elements}"
            )

            print(
                f"Selected molecular ions: "
                f"{len(selected_options)}"
            )

            print(
                f"Generated {len(candidates):,} "
                "candidate guide lines"
            )

            fig = make_spectrum_plot(
                roi=roi,
                candidates=candidates,
                lower_mass=lower_widget.value,
                upper_mass=upper_widget.value,
                bin_width=bin_width_widget.value,
                logarithmic_y=log_widget.value,
                show_labels=labels_widget.value
            )

            fig.show(
                config={
                    "scrollZoom": True,
                    "displaylogo": False,
                    "responsive": True
                }
            )

            save_path = save_widget.value.strip()

            if save_path:

                fig.write_html(
                    save_path,
                    include_plotlyjs=True,
                    full_html=True
                )

                print(
                    f"Saved interactive plot to:\n"
                    f"{save_path}"
                )

        except Exception as error:

            print(
                f"ERROR:\n"
                f"{type(error).__name__}: {error}"
            )


button.on_click(
    generate_button_clicked
)


# ============================================================
# Display GUI
# ============================================================

element_rows = []

for symbol in ALL_ELEMENTS:

    element_rows.append(
        widgets.HBox(
            [
                element_widgets[symbol],
                charge_widgets[symbol]
            ]
        )
    )


display(
    widgets.VBox(
        [
            widgets.HTML(
                "<h3>APT Candidate List Maker</h3>"
            ),

            epos_widget,

            widgets.HTML(
                "<b>Select elements and maximum atomic charge</b>"
            ),

            widgets.GridBox(
                element_rows,
                layout=widgets.Layout(
                    grid_template_columns=(
                        "repeat(3, 190px)"
                    ),
                    grid_gap="4px 12px"
                )
            ),

            widgets.HBox(
                [
                    lower_widget,
                    upper_widget,
                    bin_width_widget
                ]
            ),

            widgets.HBox(
                [
                    atomic_widget,
                    cluster_size_widget,
                    log_widget,
                    labels_widget
                ]
            ),

            widgets.HTML(
                "<b>Select individual molecular ions</b>"
            ),

            widgets.VBox(
                molecular_boxes
            ),

            manual_widget,
            save_widget,
            candidate_save_widget,
            button,
            output
        ]
    )
)
