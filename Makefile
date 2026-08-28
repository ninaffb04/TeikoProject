.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python scripts/load_data.py
	python scripts/analysis.py

dashboard:
	streamlit run scripts/dashboard.py
