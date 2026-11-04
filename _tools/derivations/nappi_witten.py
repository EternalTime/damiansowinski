"""The closed cosmology of Nappi and Witten, hep-th/9206078 (16), (21).

The metric is the large-level string-frame geometry, with its common
length scale restored.  The 1993 plane wave is a separate spacetime.
"""

METRIC = "nappi_witten"
CHARTS = ("cosmological", "bang_regular", "invariant", "whisker")
CLOCKS = {"cosmological": r"\psi", "bang_regular": "x",
          "invariant": "r", "whisker": r"\chi"}
DIMENSIONS = {(METRIC, system): {CLOCKS[system]: "1", "s": "1", r"\rho": "1",
                               (r"\lambda'" if system == "bang_regular" else r"\lambda"): "1",
                               "L": "L", "b": "1", "D": "1"}
              for system in CHARTS}


def chart_spec(system):
    """Four coordinate forms of the same classical coset geometry."""
    angular = r"\lambda'" if system == "bang_regular" else r"\lambda"
    coords = [CLOCKS[system], "s", r"\rho", angular]
    common = [r"s \in (0,\pi/2)", r"\rho \in [0,2\pi)",
              angular+r" \in [0,2\pi)", "L>0", "b>0"]
    if system == "cosmological":
        denominator = r"\sin^2s\cos^2\psi+b^2\cos^2s\sin^2\psi"
        line = (r"ds^2=L^2\left(-d\psi^2+ds^2+\dfrac{\cos^2s\cos^2\psi}{D}d\rho^2"
                r"+\dfrac{b^2\sin^2s\sin^2\psi}{D}d\lambda^2\right)")
        domains, name = [r"\psi \in (0,\pi/2)"]+common, "Cosmological"
    elif system == "bang_regular":
        denominator = r"\sin^2s+b^2x\cos^2s"
        line = (r"ds^2=L^2\left(ds^2+\dfrac{(2+x)\sin^2s-b^2\cos^2s}{4(1+x)^2D}dx^2"
                r"+\dfrac{\cos^2s}{D}d\rho^2+\dfrac{b^2x\sin^2s}{D}d\lambda'^2"
                r"+\dfrac{b\sin^2s}{D}dx\,d\lambda'\right)")
        domains, name = [r"x \in (-1,\infty)", r"D\neq0"]+common, "Across the Big Bang"
    elif system == "invariant":
        denominator = r"(1-r)\sin^2s+b^2r\cos^2s"
        line = (r"ds^2=L^2\left(-\dfrac{dr^2}{4r(1-r)}+ds^2"
                r"+\dfrac{(1-r)\cos^2s}{D}d\rho^2+\dfrac{b^2r\sin^2s}{D}d\lambda^2\right)")
        domains, name = [r"r \in (-\infty,\infty)", r"r\neq0,1", r"D\neq0"]+common, "Invariant Coordinate"
    elif system == "whisker":
        denominator = r"\sin^2s\cosh^2\chi-b^2\cos^2s\sinh^2\chi"
        line = (r"ds^2=L^2\left(d\chi^2+ds^2+\dfrac{\cos^2s\cosh^2\chi}{D}d\rho^2"
                r"-\dfrac{b^2\sin^2s\sinh^2\chi}{D}d\lambda^2\right)")
        domains, name = [r"\chi \in (0,\infty)", r"D\neq0"]+common, "Static Extension"
    else:
        raise ValueError(system)
    parameters = ["L", "b", "D = "+denominator]
    return {"metric_id": METRIC, "system": {"id": system, "name": name,
            "coords": coords, "domains": domains, "parameters": parameters,
            "line_element": line}, "chart_line_element": line, "printer": {}}
