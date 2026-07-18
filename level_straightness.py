#!/usr/bin/env python3
"""
Precision Level Straightness Workbench
======================================

Reduce precision-level step ("differential levelling") runs into a straightness
profile with live visualisation.

Method
------
1. Read the local slope of each contiguous step with the level (in divisions).
2. Cumulatively sum the readings to build the height profile h_i = sum(r).
3. Fit a reference line (endpoint / least-squares / minimum-zone).
4. Deviation e_i = h_i - line_i.  Straightness = max(e) - min(e), peak-to-valley.
5. Convert divisions to length with k = sensitivity(mm/m) * step(m).

Station 0 is the datum (reading 0); the reading at station i is the slope of the
segment ending at station i.  Station spacing is locked to the level's foot spacing
(the step length) so the cumulative sum stays honest.

Conventions follow the accompanying guide.  No silent fallbacks: every bad input is
reported with the offending row and the fix.

Dependencies: Python 3.8+, numpy, matplotlib (uses the TkAgg backend + tkinter).
Run:  python level_straightness.py
"""

import sys
import csv
import math
from dataclasses import dataclass

# ----------------------------------------------------------------------------- #
#  Guarded imports (name the cause and the fix, per house style)                 #
# ----------------------------------------------------------------------------- #
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
except Exception as exc:  # pragma: no cover
    sys.stderr.write(
        "FATAL: could not import tkinter (%s).\n"
        "Fix: install a Python build with Tk support "
        "(Debian/Ubuntu: 'sudo apt install python3-tk').\n" % exc
    )
    raise

try:
    import numpy as np
except Exception as exc:  # pragma: no cover
    sys.stderr.write(
        "FATAL: could not import numpy (%s).\nFix: 'pip install numpy'.\n" % exc
    )
    raise

try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg,
        NavigationToolbar2Tk,
    )
except Exception as exc:  # pragma: no cover
    sys.stderr.write(
        "FATAL: could not import matplotlib with the TkAgg backend (%s).\n"
        "Fix: 'pip install matplotlib'.\n" % exc
    )
    raise


ARCSEC_PER_MM_PER_M = 206.2648  # 1 mm/m of slope = 206.2648 arc-seconds

METHODS = ("endpoint", "least_squares", "minimum_zone")
METHOD_LABELS = {
    "endpoint": "Endpoint line",
    "least_squares": "Least-squares mean line",
    "minimum_zone": "Minimum zone (ISO 12780)",
}


# ============================================================================= #
#  COMPUTATION CORE  (no GUI; independently testable)                            #
# ============================================================================= #
@dataclass
class LevelConfig:
    """Instrument + geometry that turns divisions into length."""
    sensitivity_mm_per_m: float   # slope represented by one division
    step_length_mm: float         # foot spacing == station spacing

    @property
    def k_mm_per_div(self) -> float:
        """Rise in mm produced by one division over one step."""
        return self.sensitivity_mm_per_m * (self.step_length_mm / 1000.0)

    @property
    def k_um_per_div(self) -> float:
        return self.k_mm_per_div * 1000.0


def endpoint_line(x, h):
    """Reference line through the first and last profile points -> (slope, intercept)."""
    if x[-1] == x[0]:
        raise ValueError("Endpoint line undefined: first and last station coincide.")
    m = (h[-1] - h[0]) / (x[-1] - x[0])
    b = h[0] - m * x[0]
    return m, b


def least_squares_line(x, h):
    """Best-fit line minimising sum of squared deviations -> (slope, intercept)."""
    m, b = np.polyfit(x, h, 1)
    return float(m), float(b)


def minimum_zone_line(x, h, iterations=200):
    """
    Minimum-zone reference: the slope that minimises the peak-to-valley of the
    residuals.  P-V(slope) is convex, so a ternary search converges reliably.
    Returns (slope, intercept) with the intercept placed at the zone mid-line.
    """
    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)

    def pv(m):
        r = h - m * x
        return r.max() - r.min()

    # The optimum slope lies within the span of local segment slopes; bracket wide.
    seg_slopes = np.diff(h) / np.diff(x)
    lo = float(seg_slopes.min())
    hi = float(seg_slopes.max())
    if lo == hi:                       # perfectly straight run
        lo -= 1e-9
        hi += 1e-9

    for _ in range(iterations):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        if pv(m1) < pv(m2):
            hi = m2
        else:
            lo = m1
    m = 0.5 * (lo + hi)
    r = h - m * x
    b = 0.5 * (r.max() + r.min())      # centre the zone on the line
    return float(m), float(b)


def classify_shape(dev, tol_div):
    """Describe the deviation profile: concave / convex / wavy / straight."""
    span = float(dev.max() - dev.min())
    if span <= tol_div:
        return "straight (within resolution)"
    signs = np.sign(np.where(np.abs(dev) < tol_div, 0.0, dev))
    signs = signs[signs != 0]
    if signs.size == 0:
        return "straight (within resolution)"
    changes = int(np.count_nonzero(np.diff(signs) != 0))
    mean = float(dev.mean())
    if changes <= 1:
        return "concave / hollow" if mean < 0 else "convex / crowned"
    return "wavy (%d undulations)" % (changes + 1)


@dataclass
class Reduction:
    """Everything the plot and the report need."""
    x: np.ndarray
    readings: np.ndarray
    cumulative: np.ndarray
    ref_line: np.ndarray
    deviation: np.ndarray
    slope_div_per_mm: float
    intercept_div: float
    method: str
    config: LevelConfig

    @property
    def straightness_div(self):
        return float(self.deviation.max() - self.deviation.min())

    @property
    def straightness_um(self):
        return self.straightness_div * self.config.k_um_per_div

    @property
    def dev_um(self):
        return self.deviation * self.config.k_um_per_div

    @property
    def rms_um(self):
        return float(np.sqrt(np.mean(self.dev_um ** 2)))

    @property
    def overall_tilt_mm_per_m(self):
        # end-to-end inclination of the surface in mm/m
        length_m = (self.x[-1] - self.x[0]) / 1000.0
        if length_m == 0:
            return 0.0
        rise_mm = (self.cumulative[-1] - self.cumulative[0]) * self.config.k_mm_per_div
        return rise_mm / length_m

    def shape(self, tol_frac=0.05):
        # Describe the surface shape against the endpoint chord so the label is
        # method-independent (relative to a mean line, any real profile crosses it
        # and would read as "wavy").  tol_frac ignores crossings below that
        # fraction of the P-V span, so a smooth hollow is not called wavy.
        m, b = endpoint_line(self.x, self.cumulative)
        chord_dev = self.cumulative - (m * self.x + b)
        span = float(chord_dev.max() - chord_dev.min())
        tol = max(1e-9, tol_frac * span)
        return classify_shape(chord_dev, tol)


def reduce_run(readings, config: LevelConfig, method: str) -> Reduction:
    """Full reduction of one run.  readings[0] is the datum (station 0)."""
    if method not in METHODS:
        raise ValueError("Unknown method %r; expected one of %s." % (method, METHODS))
    readings = np.asarray(readings, dtype=float)
    n = readings.size
    if n < 2:
        raise ValueError("Need at least 2 stations (datum + 1 reading); got %d." % n)

    x = np.arange(n, dtype=float) * config.step_length_mm
    cumulative = np.cumsum(readings)

    if method == "endpoint":
        m, b = endpoint_line(x, cumulative)
    elif method == "least_squares":
        m, b = least_squares_line(x, cumulative)
    else:
        m, b = minimum_zone_line(x, cumulative)

    ref = m * x + b
    dev = cumulative - ref
    return Reduction(
        x=x, readings=readings, cumulative=cumulative, ref_line=ref,
        deviation=dev, slope_div_per_mm=m, intercept_div=b,
        method=method, config=config,
    )


# ============================================================================= #
#  GUI                                                                           #
# ============================================================================= #
EXAMPLE_READINGS = [0.0, 3.5, 4.0, 5.5, 4.0, 3.0, 6.0, 5.0]  # from the guide


class StraightnessApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Precision Level Straightness Workbench")
        self.geometry("1180x760")
        self.minsize(980, 640)

        self.row_vars = []          # list of (station_var, reading_var, reversed_var)
        self.use_reversal = tk.BooleanVar(value=False)
        self.method_var = tk.StringVar(value="least_squares")
        self.sens_var = tk.StringVar(value="0.02")
        self.sens_unit_var = tk.StringVar(value="mm/m")
        self.step_var = tk.StringVar(value="100")
        self.nseg_var = tk.IntVar(value=len(EXAMPLE_READINGS) - 1)
        self.held_traces = []       # snapshots for overlay comparison
        self._last_reduction = None

        self._build_layout()
        self._rebuild_table(initial=EXAMPLE_READINGS)
        self.log("Loaded example run from the guide. Press Update to reduce.")
        self.update_reduction()

    # ---------------------------------------------------------------- layout -- #
    def _build_layout(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # ---- instrument config ------------------------------------------------
        cfg = ttk.LabelFrame(left, text="Instrument", padding=8)
        cfg.pack(fill=tk.X)
        ttk.Label(cfg, text="Sensitivity / division").grid(row=0, column=0, sticky="w")
        ttk.Entry(cfg, textvariable=self.sens_var, width=9).grid(row=0, column=1, padx=4)
        unit = ttk.Combobox(cfg, textvariable=self.sens_unit_var, width=8,
                            state="readonly", values=("mm/m", "arc-sec"))
        unit.grid(row=0, column=2)
        unit.bind("<<ComboboxSelected>>", lambda e: self.update_reduction())
        ttk.Label(cfg, text="Step / foot spacing (mm)").grid(row=1, column=0, sticky="w",
                                                            pady=(6, 0))
        ttk.Entry(cfg, textvariable=self.step_var, width=9).grid(row=1, column=1,
                                                                padx=4, pady=(6, 0))
        ttk.Label(cfg, text="Segments").grid(row=2, column=0, sticky="w", pady=(6, 0))
        seg = ttk.Spinbox(cfg, from_=1, to=200, textvariable=self.nseg_var, width=7,
                        command=self._on_segments_changed)
        seg.grid(row=2, column=1, padx=4, pady=(6, 0), sticky="w")

        # ---- reference method -------------------------------------------------
        meth = ttk.LabelFrame(left, text="Reference line", padding=8)
        meth.pack(fill=tk.X, pady=(8, 0))
        for m in METHODS:
            ttk.Radiobutton(meth, text=METHOD_LABELS[m], value=m,
                            variable=self.method_var,
                            command=self.update_reduction).pack(anchor="w")
        ttk.Checkbutton(meth, text="Use reversal (2-reading zero cancel)",
                        variable=self.use_reversal,
                        command=self._toggle_reversal).pack(anchor="w", pady=(6, 0))

        # ---- readings table ---------------------------------------------------
        tbl = ttk.LabelFrame(left, text="Readings (divisions)", padding=6)
        tbl.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        head = ttk.Frame(tbl)
        head.pack(fill=tk.X)
        ttk.Label(head, text="Station (mm)", width=12,
                anchor="center").grid(row=0, column=0)
        ttk.Label(head, text="Reading", width=10, anchor="center").grid(row=0, column=1)
        self.rev_head = ttk.Label(head, text="Reversed", width=10, anchor="center")
        # reversed header shown only when reversal enabled

        canvas = tk.Canvas(tbl, highlightthickness=0, width=280)
        scroll = ttk.Scrollbar(tbl, orient="vertical", command=canvas.yview)
        self.rows_frame = ttk.Frame(canvas)
        self.rows_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.rows_frame, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # ---- action buttons ---------------------------------------------------
        btns = ttk.Frame(left)
        btns.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(btns, text="Update", command=self.update_reduction).grid(
            row=0, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Hold trace", command=self.hold_trace).grid(
            row=0, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Clear traces", command=self.clear_traces).grid(
            row=0, column=2, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Load CSV", command=self.load_csv).grid(
            row=1, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Save CSV", command=self.save_csv).grid(
            row=1, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Load example", command=self.load_example).grid(
            row=1, column=2, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Save plot (PNG)", command=self.save_plot).grid(
            row=2, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(btns, text="Export report", command=self.export_report).grid(
            row=2, column=1, columnspan=2, sticky="ew", padx=2, pady=2)
        for c in range(3):
            btns.columnconfigure(c, weight=1)

        # ---- results ----------------------------------------------------------
        res = ttk.LabelFrame(right, text="Result", padding=6)
        res.pack(fill=tk.X)
        self.result_text = tk.Text(res, height=6, wrap="word", font=("TkFixedFont", 10))
        self.result_text.pack(fill=tk.X)
        self.result_text.configure(state="disabled")

        # ---- figure -----------------------------------------------------------
        self.fig = Figure(figsize=(7.4, 5.4), dpi=100, constrained_layout=True)
        self.ax_profile = self.fig.add_subplot(2, 1, 1)
        self.ax_dev = self.fig.add_subplot(2, 1, 2)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        NavigationToolbar2Tk(self.canvas, right).update()

        # ---- status / log -----------------------------------------------------
        self.status = tk.StringVar(value="Ready.")
        bar = ttk.Frame(self)
        bar.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(bar, textvariable=self.status, relief="sunken",
                anchor="w", padding=3).pack(fill=tk.X)

    # ------------------------------------------------------------- table mgmt -- #
    def _make_row(self, station, reading, reversed_val=""):
        sv = tk.StringVar(value=("%g" % station))
        rv = tk.StringVar(value=("" if reading is None else "%g" % reading))
        xv = tk.StringVar(value=("" if reversed_val in ("", None) else "%g" % reversed_val))
        idx = len(self.row_vars)
        r = ttk.Frame(self.rows_frame)
        r.grid(row=idx, column=0, sticky="w")
        ttk.Label(r, textvariable=sv, width=12, anchor="center",
                relief="groove").grid(row=0, column=0, padx=1, pady=1)
        e_read = ttk.Entry(r, textvariable=rv, width=10, justify="center")
        e_read.grid(row=0, column=1, padx=1, pady=1)
        e_read.bind("<Return>", lambda ev: self.update_reduction())
        e_read.bind("<FocusOut>", lambda ev: self.update_reduction())
        e_rev = ttk.Entry(r, textvariable=xv, width=10, justify="center")
        e_rev.bind("<Return>", lambda ev: self.update_reduction())
        e_rev.bind("<FocusOut>", lambda ev: self.update_reduction())
        if self.use_reversal.get():
            e_rev.grid(row=0, column=2, padx=1, pady=1)
        self.row_vars.append((sv, rv, xv))
        self._row_widgets = getattr(self, "_row_widgets", [])
        self._row_widgets.append((r, e_rev))

    def _rebuild_table(self, initial=None):
        for child in self.rows_frame.winfo_children():
            child.destroy()
        self.row_vars = []
        self._row_widgets = []
        try:
            step = float(self.step_var.get())
        except ValueError:
            step = 100.0
        nseg = max(1, int(self.nseg_var.get()))
        for i in range(nseg + 1):
            reading = None
            if initial is not None and i < len(initial):
                reading = initial[i]
            elif i == 0:
                reading = 0.0
            self._make_row(i * step, reading)
        # station 0 is the datum
        if self.row_vars:
            self.row_vars[0][1].set("0")

    def _on_segments_changed(self):
        existing = [rv.get() for (_, rv, _) in self.row_vars]
        vals = []
        for v in existing:
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(None)
        self._rebuild_table(initial=vals)
        self.update_reduction()

    def _toggle_reversal(self):
        # show/hide the reversed header + entry column
        if self.use_reversal.get():
            self.rev_head.grid(row=0, column=2)
        else:
            self.rev_head.grid_forget()
        for (frame, e_rev) in getattr(self, "_row_widgets", []):
            if self.use_reversal.get():
                e_rev.grid(row=0, column=2, padx=1, pady=1)
            else:
                e_rev.grid_forget()
        self.update_reduction()

    # -------------------------------------------------------------- reducing -- #
    def _current_config(self):
        try:
            sens = float(self.sens_var.get())
        except ValueError:
            raise ValueError("Sensitivity %r is not a number." % self.sens_var.get())
        if sens <= 0:
            raise ValueError("Sensitivity must be > 0; got %g." % sens)
        if self.sens_unit_var.get() == "arc-sec":
            sens = sens / ARCSEC_PER_MM_PER_M   # convert to mm/m
        try:
            step = float(self.step_var.get())
        except ValueError:
            raise ValueError("Step length %r is not a number." % self.step_var.get())
        if step <= 0:
            raise ValueError("Step length must be > 0; got %g." % step)
        return LevelConfig(sensitivity_mm_per_m=sens, step_length_mm=step)

    def _gather_readings(self):
        readings = []
        for i, (sv, rv, xv) in enumerate(self.row_vars):
            raw = rv.get().strip()
            if raw == "":
                if i == 0:
                    readings.append(0.0)
                    continue
                raise ValueError(
                    "Station %s mm (row %d) reading is blank. "
                    "Enter a number, or reduce the segment count." % (sv.get(), i))
            try:
                val = float(raw)
            except ValueError:
                raise ValueError(
                    "Station %s mm (row %d) reading %r is not a number." %
                    (sv.get(), i, raw))
            if self.use_reversal.get():
                rraw = xv.get().strip()
                if rraw != "":
                    try:
                        rev = float(rraw)
                    except ValueError:
                        raise ValueError(
                            "Station %s mm (row %d) reversed reading %r is not a "
                            "number." % (sv.get(), i, rraw))
                    val = (val - rev) / 2.0     # zero-cancelled true slope
            readings.append(val)
        return readings

    def update_reduction(self):
        try:
            config = self._current_config()
            readings = self._gather_readings()
            reduction = reduce_run(readings, config, self.method_var.get())
        except ValueError as exc:
            self.log("ERROR: %s" % exc, error=True)
            return
        self._last_reduction = reduction
        self._draw(reduction)
        self._show_result(reduction)
        self.log("Reduced %d stations, method=%s, straightness=%.3f div (%.2f um)."
                % (reduction.x.size, reduction.method,
                    reduction.straightness_div, reduction.straightness_um))

    # -------------------------------------------------------------- plotting -- #
    def _draw(self, red: Reduction):
        self.ax_profile.clear()
        self.ax_dev.clear()

        # profile + reference line
        self.ax_profile.plot(red.x, red.cumulative, "o-", color="#1f77b4",
                            label="measured profile")
        self.ax_profile.plot(red.x, red.ref_line, "--", color="#d62728",
                            label="%s" % METHOD_LABELS[red.method])
        self.ax_profile.set_ylabel("cumulative height (div)")
        self.ax_profile.set_title("Height profile vs reference line")
        self.ax_profile.legend(fontsize=8, loc="best")
        self.ax_profile.grid(True, alpha=0.3)

        # deviation (straightness) with zone band, in micrometres
        dev_um = red.dev_um
        self.ax_dev.axhspan(dev_um.min(), dev_um.max(), color="#d62728", alpha=0.08)
        self.ax_dev.axhline(dev_um.max(), color="#d62728", lw=0.8, ls=":")
        self.ax_dev.axhline(dev_um.min(), color="#d62728", lw=0.8, ls=":")
        for trace in self.held_traces:
            self.ax_dev.plot(trace["x"], trace["dev_um"], "-",
                            color="0.7", lw=1.0, alpha=0.8)
        self.ax_dev.plot(red.x, dev_um, "o-", color="#2ca02c",
                        label="deviation (P-V = %.2f um)" % red.straightness_um)
        self.ax_dev.axhline(0, color="k", lw=0.6)
        # mark extremes
        i_max = int(np.argmax(dev_um))
        i_min = int(np.argmin(dev_um))
        self.ax_dev.annotate("max", (red.x[i_max], dev_um[i_max]),
                            textcoords="offset points", xytext=(0, 6), fontsize=8)
        self.ax_dev.annotate("min", (red.x[i_min], dev_um[i_min]),
                            textcoords="offset points", xytext=(0, -12), fontsize=8)
        self.ax_dev.set_xlabel("station (mm)")
        self.ax_dev.set_ylabel("deviation (um)")
        self.ax_dev.set_title("Straightness deviation")
        self.ax_dev.legend(fontsize=8, loc="best")
        self.ax_dev.grid(True, alpha=0.3)
        self.canvas.draw_idle()

    # --------------------------------------------------------------- results -- #
    def _show_result(self, red: Reduction):
        cfg = red.config
        lines = [
            "Straightness (peak-to-valley) : %8.3f div   = %8.2f um" % (
                red.straightness_div, red.straightness_um),
            "RMS deviation                 : %8.2f um" % red.rms_um,
            "Reference slope               : %10.5f div/mm" % red.slope_div_per_mm,
            "Overall tilt                  : %8.3f mm/m" % red.overall_tilt_mm_per_m,
            "Shape                         : %s" % red.shape(),
            "k (rise per division per step): %8.3f um    [%s, step %g mm]" % (
                cfg.k_um_per_div, "%.4g mm/m" % cfg.sensitivity_mm_per_m,
                cfg.step_length_mm),
        ]
        self.result_text.configure(state="normal")
        self.result_text.delete("1.0", tk.END)
        self.result_text.insert("1.0", "\n".join(lines))
        self.result_text.configure(state="disabled")

    # ----------------------------------------------------------------- traces -- #
    def hold_trace(self):
        if self._last_reduction is None:
            self.log("Nothing to hold: reduce a run first.", error=True)
            return
        red = self._last_reduction
        self.held_traces.append({"x": red.x.copy(), "dev_um": red.dev_um.copy()})
        self.log("Held trace #%d for comparison." % len(self.held_traces))
        self._draw(red)

    def clear_traces(self):
        self.held_traces = []
        self.log("Cleared held traces.")
        if self._last_reduction is not None:
            self._draw(self._last_reduction)

    # -------------------------------------------------------------------- I/O -- #
    def load_example(self):
        self.step_var.set("100")
        self.nseg_var.set(len(EXAMPLE_READINGS) - 1)
        self._rebuild_table(initial=EXAMPLE_READINGS)
        self.log("Loaded example run.")
        self.update_reduction()

    def load_csv(self):
        path = filedialog.askopenfilename(
            title="Load readings CSV",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        try:
            stations, readings = [], []
            with open(path, newline="") as fh:
                reader = csv.reader(fh)
                for row in reader:
                    if not row or not row[0].strip():
                        continue
                    cell = row[0].strip().lower()
                    if cell in ("station", "x", "pos", "position"):
                        continue                      # header
                    stations.append(float(row[0]))
                    readings.append(float(row[1]) if len(row) > 1 and
                                    row[1].strip() != "" else 0.0)
            if len(readings) < 2:
                raise ValueError("Need at least 2 rows; found %d." % len(readings))
            step = stations[1] - stations[0] if len(stations) > 1 else 100.0
            if step <= 0:
                raise ValueError("Station spacing must increase; got step %g." % step)
            self.step_var.set("%g" % step)
            self.nseg_var.set(len(readings) - 1)
            self._rebuild_table(initial=readings)
            self.update_reduction()
            self.log("Loaded %d rows from %s (step %g mm)." %
                    (len(readings), path, step))
        except Exception as exc:
            self.log("ERROR loading CSV: %s" % exc, error=True)
            messagebox.showerror("Load CSV", str(exc))

    def save_csv(self):
        if self._last_reduction is None:
            self.log("Nothing to save: reduce a run first.", error=True)
            return
        path = filedialog.asksaveasfilename(
            title="Save reduction CSV", defaultextension=".csv",
            filetypes=[("CSV", "*.csv")])
        if not path:
            return
        red = self._last_reduction
        try:
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["station_mm", "reading_div", "cumulative_div",
                            "reference_div", "deviation_div", "deviation_um"])
                for i in range(red.x.size):
                    w.writerow(["%g" % red.x[i], "%g" % red.readings[i],
                                "%g" % red.cumulative[i], "%.6g" % red.ref_line[i],
                                "%.6g" % red.deviation[i], "%.6g" % red.dev_um[i]])
            self.log("Saved reduction to %s." % path)
        except Exception as exc:
            self.log("ERROR saving CSV: %s" % exc, error=True)
            messagebox.showerror("Save CSV", str(exc))

    def save_plot(self):
        path = filedialog.asksaveasfilename(
            title="Save plot", defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf")])
        if not path:
            return
        try:
            self.fig.savefig(path, dpi=150)
            self.log("Saved plot to %s." % path)
        except Exception as exc:
            self.log("ERROR saving plot: %s" % exc, error=True)

    def export_report(self):
        if self._last_reduction is None:
            self.log("Nothing to export: reduce a run first.", error=True)
            return
        path = filedialog.asksaveasfilename(
            title="Export report", defaultextension=".txt",
            filetypes=[("Text", "*.txt")])
        if not path:
            return
        red = self._last_reduction
        cfg = red.config
        try:
            with open(path, "w") as fh:
                fh.write("Precision Level Straightness Report\n")
                fh.write("===================================\n\n")
                fh.write("Instrument sensitivity : %.4g mm/m/div\n"
                        % cfg.sensitivity_mm_per_m)
                fh.write("Step / foot spacing    : %g mm\n" % cfg.step_length_mm)
                fh.write("k (um per div per step): %.3f um\n" % cfg.k_um_per_div)
                fh.write("Reference method       : %s\n" % METHOD_LABELS[red.method])
                fh.write("Reference slope        : %.6g div/mm\n"
                        % red.slope_div_per_mm)
                fh.write("\nStraightness (P-V)     : %.3f div = %.2f um\n"
                        % (red.straightness_div, red.straightness_um))
                fh.write("RMS deviation          : %.2f um\n" % red.rms_um)
                fh.write("Overall tilt           : %.3f mm/m\n"
                        % red.overall_tilt_mm_per_m)
                fh.write("Shape                  : %s\n\n" % red.shape())
                fh.write("%10s %10s %12s %12s %12s\n" %
                        ("stn(mm)", "read", "cumul", "dev(div)", "dev(um)"))
                for i in range(red.x.size):
                    fh.write("%10g %10g %12.4f %12.4f %12.3f\n" %
                            (red.x[i], red.readings[i], red.cumulative[i],
                            red.deviation[i], red.dev_um[i]))
            self.log("Exported report to %s." % path)
        except Exception as exc:
            self.log("ERROR exporting report: %s" % exc, error=True)

    # ----------------------------------------------------------------- status -- #
    def log(self, msg, error=False):
        self.status.set(msg)
        stream = sys.stderr if error else sys.stdout
        stream.write(msg + "\n")


def main():
    app = StraightnessApp()
    app.mainloop()


if __name__ == "__main__":
    main()
