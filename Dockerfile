# Dockerfile for the GenoDynLabs backend on Railway.
#
# Ships the prebuilt data/reference.pkl (a few MB) so the image never parses the
# 33 MB SNPedia CSVs at build or run time — smaller image, lower memory, faster
# cold start. Regenerate the pkl locally with the Makefile when the data changes.
#
# The SAME image runs both services; Railway just overrides the start command:
#   web    -> gunicorn app:app ...     (default CMD below)
#   worker -> python worker.py         (set as the worker service's start cmd)

FROM python:3.11-slim-bookworm

# System libraries WeasyPrint needs (Pango/Cairo/etc.) plus fonts.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpango-1.0-0 \
    libpangoft2-1.0-0 \
    libharfbuzz0b \
    libffi-dev \
    libcairo2 \
    libgdk-pixbuf2.0-0 \
    libglib2.0-0 \
    shared-mime-info \
    fonts-liberation \
    fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PYTHONUNBUFFERED=1
EXPOSE 8080

# Default = web. Worker service overrides this with: python worker.py
CMD ["sh", "-c", "gunicorn app:app --bind 0.0.0.0:${PORT:-8080} --timeout 120 --workers ${WEB_WORKERS:-2} --threads ${WEB_THREADS:-4}"]
