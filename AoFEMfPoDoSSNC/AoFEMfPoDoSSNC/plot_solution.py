import ast
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT_FILE = ROOT / "output_for_plot.txt"
OUTPUT_IMAGE = ROOT / "solution_plot.png"


ALLOWED_FUNCTIONS = {
    "sin": np.sin,
    "cos": np.cos,
    "tg": np.tan,
    "tan": np.tan,
    "ctg": lambda x: 1.0 / np.tan(x),
    "asin": np.arcsin,
    "acos": np.arccos,
    "atg": np.arctan,
    "atan": np.arctan,
    "actg": lambda x: np.pi / 2.0 - np.arctan(x),
    "sqrt": np.sqrt,
    "sqr": lambda x: x * x,
    "exp": np.exp,
    "log": np.log,
    "abs": np.abs,
}

ALLOWED_CONSTANTS = {
    "pi": np.pi,
    "e": math.e,
}


def piecewise(*args):
    if len(args) < 3 or len(args) % 2 == 0:
        raise ValueError("pw expects: pw(border1, expr1, ..., lastExpr)")

    result = np.array(args[-1], copy=True)
    for index in range(len(args) - 3, -1, -2):
        border = args[index]
        expr = args[index + 1]
        result = np.where(piecewise._x <= border, expr, result)
    return result


ALLOWED_FUNCTIONS["piecewise"] = piecewise
ALLOWED_FUNCTIONS["pw"] = piecewise


def read_plot_data(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")

    data = {
        "meta": {},
        "coefficients": [],
        "pieces": [],
        "nodes": [],
    }
    section = "meta"

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue

        if section == "meta":
            key, value = line.split(" ", 1)
            data["meta"][key] = value
        elif section == "coefficients":
            index, value = line.split()
            data["coefficients"].append((int(index), float(value)))
        elif section == "pieces":
            left, right, free_term, x_coeff = map(float, line.split())
            data["pieces"].append((left, right, free_term, x_coeff))
        elif section == "nodes":
            x, value = map(float, line.split())
            data["nodes"].append((x, value))

    return data


def validate_expression(node):
    allowed_nodes = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Call,
        ast.Name,
        ast.Load,
        ast.Constant,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.Pow,
        ast.USub,
        ast.UAdd,
    )
    if not isinstance(node, allowed_nodes):
        raise ValueError(f"Unsupported expression node: {type(node).__name__}")

    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in ALLOWED_FUNCTIONS:
            raise ValueError("Unsupported function in exact_u")

    if isinstance(node, ast.Name):
        if node.id != "x" and node.id not in ALLOWED_CONSTANTS and node.id not in ALLOWED_FUNCTIONS:
            raise ValueError(f"Unsupported name in exact_u: {node.id}")

    for child in ast.iter_child_nodes(node):
        validate_expression(child)


def compile_expression(expression):
    expression = expression.replace("^", "**")
    tree = ast.parse(expression, mode="eval")
    validate_expression(tree)
    code = compile(tree, "<exact_u>", "eval")

    def evaluator(x):
        piecewise._x = x
        namespace = {"x": x, **ALLOWED_FUNCTIONS, **ALLOWED_CONSTANTS}
        return eval(code, {"__builtins__": {}}, namespace)

    return evaluator


def sample_piecewise_solution(pieces, points_per_piece=40):
    xs = []
    ys = []
    for left, right, free_term, x_coeff in pieces:
        local_x = np.linspace(left, right, points_per_piece, endpoint=True)
        local_y = free_term + x_coeff * local_x
        xs.append(local_x)
        ys.append(local_y)

    return np.concatenate(xs), np.concatenate(ys)


def evaluate_piecewise_solution(pieces, x):
    result = np.zeros_like(x, dtype=float)
    for index, (left, right, free_term, x_coeff) in enumerate(pieces):
        if index == len(pieces) - 1:
            mask = (x >= left) & (x <= right)
        else:
            mask = (x >= left) & (x < right)
        result[mask] = free_term + x_coeff * x[mask]
    return result


def main():
    data = read_plot_data(INPUT_FILE)
    meta = data["meta"]
    pieces = data["pieces"]
    nodes = data["nodes"]

    if not pieces or not nodes:
        raise ValueError("Plot file does not contain pieces/nodes")

    x_approx, y_approx = sample_piecewise_solution(pieces)
    x_nodes = np.array([item[0] for item in nodes])
    y_nodes = np.array([item[1] for item in nodes])

    fig, ax = plt.subplots(figsize=(10, 6), dpi=140)
    ax.plot(x_approx, y_approx, color="#1f77b4", linewidth=2.2, label="u_h(x), FEM")
    ax.scatter(x_nodes, y_nodes, color="#1f77b4", s=24, zorder=3, label="mesh nodes")

    exact_expression = meta.get("exact_u", "").strip()
    if exact_expression:
        exact = compile_expression(exact_expression)
        x_exact = np.linspace(float(meta["l"]) * 0.0, float(meta["l"]), 800)
        y_exact = exact(x_exact)
        y_error = evaluate_piecewise_solution(pieces, x_exact) - y_exact
        max_error = float(np.max(np.abs(y_error)))
        exact_label = "exact u(x)"
        if len(exact_expression) <= 45:
            exact_label = f"{exact_label} = {exact_expression}"
        ax.plot(x_exact, y_exact, color="#d62728", linewidth=1.8, linestyle="--", label=exact_label)
        print(f"Max |u_h(x)-u_exact(x)| on plot grid: {max_error:.6e}")

    limiter = float(meta["m"])
    ax.axhline(limiter, color="#555555", linewidth=1.0, linestyle=":", label="limiter")
    ax.axhline(-limiter, color="#555555", linewidth=1.0, linestyle=":")

    title_prefix = meta.get("name", "").strip()
    title = f"N={meta.get('N')}, l={meta.get('l')}, m={meta.get('m')}"
    if title_prefix:
        title = f"{title_prefix}\n{title}"
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("u(x)")
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_IMAGE)
    plt.close(fig)

    print(f"Saved plot to {OUTPUT_IMAGE}")


if __name__ == "__main__":
    main()
