# Reproducible commands. `make check` runs everything CI runs.
.PHONY: install test gate lint check schema legend rd r-check

install:
	pip install -e ".[dev]"

test:
	pytest -q

gate:
	python checks/gate.py

lint:
	ruff check .

check: lint gate test

# the R companion package (r/): built and checked the way CRAN does, where R is installed
r-check:
	R CMD build r
	R CMD check --no-manual foodnet_*.tar.gz
	rm -rf foodnet_*.tar.gz foodnet.Rcheck

schema:
	python -m foodnet schema --out schema/metabolite_network.schema.json

legend:
	python -c "from foodnet.legend import legend_svg; open('docs/legend.svg', 'w').write(legend_svg())"
	python -c "from foodnet.brand import logo_svg; open('docs/logo.svg', 'w').write(logo_svg(64))"

# the R package's help pages, from the #' comments in r/R (no roxygen2 needed)
rd:
	python packaging/rd.py
