# Dockerfile for GenoDynLabs backend on Railway
#
# WeasyPrint needs system libraries (Pango, Cairo, GObject, etc.) that aren't
# present in a default Python image, so we install them explicitly here.

FROM python:3.11-slim

# System libraries required by WeasyPrint + general image build
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
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

# Install Python deps first (better Docker layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Railway provides $PORT at runtime
ENV PYTHONUNBUFFERED=1
EXPOSE 8080

CMD gunicorn app:app --bind 0.0.0.0:${PORT:-8080} --timeout 300 --workers 2
