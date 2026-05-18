# Computational Mathematics
## :construction: Work In Progress :construction:
My attempts at implementing algorithms relating to computational mathematics \
provided as a python package.

Key features:
- :zap: probably not very fast

- :white_check_mark: probably not very correct

- :crab: definitely handwritten in rust

<p>¯\_(ツ)_/¯</p>

## Setup
Current set up for this project requires both python and rust and uses the tools `uv` and `maturin`. \

Steps:
- Install uv: `pip install uv`

- Install maturin: `uv tool install maturin`

- Create venv: `uv sync`

- Start developing with maturin: `maturin develop` or `maturin develop -r`

There are two sets of sets, in rust and python. \
These can be run with `cargo test` and `pytest` respectively.
