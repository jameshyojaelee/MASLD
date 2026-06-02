.PHONY: run build-parquet check deploy deploy-data deploy-space

# PYTHONNOUSERSITE=1 isolates from ~/.local site-packages so that a pip --user
# install of huggingface_hub / typer / numexpr does not clash with the active
# micromamba env's numpy/pandas ABI.
run:
	PYTHONNOUSERSITE=1 streamlit run app.py

build-parquet:
	python -m data.build_parquet

check:
	python -c "from data.paths import data_dir; print('DATA_DIR =', data_dir()); p = data_dir() / 'atlas_summary.json'; assert p.exists(), f'missing {p}'; print('OK')"

deploy-data:
	./deploy/upload-data.sh

deploy-space:
	./deploy/push-space.sh

deploy: deploy-data deploy-space
