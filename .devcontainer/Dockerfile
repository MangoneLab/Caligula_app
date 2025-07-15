FROM python:3.10

WORKDIR /workspace

COPY ../requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY ../caligula_single.py .
COPY ../extract_strand_features_fixed.py .
COPY ../model.pkl .
