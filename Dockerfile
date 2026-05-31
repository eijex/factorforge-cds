FROM python:3.11-slim

WORKDIR /app

# Install package dependencies
COPY pyproject.toml .
COPY src/ src/

RUN pip install --no-cache-dir -e .

# Copy runtime files
COPY api/ api/
COPY web/ web/
COPY data/ data/
COPY scripts/serve.py .

EXPOSE 8080

CMD ["python", "serve.py"]
