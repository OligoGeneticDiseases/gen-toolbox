# Docker Deployment Guide

This guide explains how to build and deploy the Oligogenicity Analysis Toolbox using Docker.

## Quick Start

### Building the Image

```bash
# Build the Docker image
docker build -t oligogenicity:latest .

# Or build with docker-compose
docker-compose build
```

### Running with Docker

```bash
# Show help
docker run --rm oligogenicity:latest

# Run a specific command with mounted volumes
docker run --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/tmp:/app/tmp \
  oligogenicity:latest findtype \
  --source /app/data \
  --directory /app/output \
  --type .vcf
```

### Running with Docker Compose

```bash
# Run the main service
docker-compose up gen-toolbox

# Run in development mode
docker-compose --profile dev up gen-toolbox-dev

# Run a one-off command
docker-compose run --rm gen-toolbox findtype --source /app/data --directory /app/output --type .vcf
```

## Volume Mounts

### Required Directories

Create these directories on your host system:

```bash
mkdir -p data output tmp
```

### Data Directory Structure

```
data/
├── vcfs/           # Input VCF files
├── metadata.txt    # Sample metadata file
└── vep_cache/      # VEP cache (if using local VEP)

output/             # Results will be written here
tmp/                # Temporary files during processing
```

## Configuration

### Environment Variables

- `JAVA_OPTS`: JVM options for Spark (default: `-Xmx8g`)
- `SPARK_LOCAL_DIRS`: Temporary directory for Spark operations

### VEP Configuration

If using VEP annotation, update the volume mount in `docker-compose.yml`:

```yaml
volumes:
  - /your/vep/cache/path:/vep/cache:ro
```

And update `src/config/vep_settings.json` accordingly.

## Example Commands

### Find VCF Files

```bash
docker run --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/output:/app/output \
  oligogenicity:latest findtype \
  --source /app/data \
  --directory /app/output \
  --type .vcf
```

### Process VCF Files

```bash
docker run --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/tmp:/app/tmp \
  oligogenicity:latest readvcfs \
  --file /app/data/vcfs \
  --dest /app/output \
  --globals /app/data/metadata.txt \
  --phenotype "Case" \
  --temp /app/tmp
```

### Run SKAT Analysis

```bash
docker run --rm \
  -v $(pwd)/data:/app/data:ro \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/tmp:/app/tmp \
  oligogenicity:latest runskat \
  --file /app/data/vcfs \
  --dest /app/output \
  --globals /app/data/metadata.txt \
  --phenotype "Case" \
  --temp /app/tmp
```

## Publishing to Registry

### Docker Hub

```bash
# Tag the image
docker tag oligogenicity:latest your-username/oligogenicity:latest
docker tag oligogenicity:latest your-username/oligogenicity:1.0.0

# Push to Docker Hub
docker push your-username/oligogenicity:latest
docker push your-username/oligogenicity:1.0.0
```

### GitHub Container Registry

```bash
# Tag for GitHub Container Registry
docker tag oligogenicity:latest ghcr.io/your-username/oligogenicity:latest
docker tag oligogenicity:latest ghcr.io/your-username/oligogenicity:1.0.0

# Login and push
echo $GITHUB_TOKEN | docker login ghcr.io -u your-username --password-stdin
docker push ghcr.io/your-username/oligogenicity:latest
docker push ghcr.io/your-username/oligogenicity:1.0.0
```

## Resource Requirements

### Minimum Requirements
- **Memory:** 8GB RAM
- **CPU:** 4 cores
- **Disk:** 50GB free space

### Recommended for Large Datasets
- **Memory:** 32GB+ RAM
- **CPU:** 16+ cores
- **Disk:** 500GB+ free space

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**
   ```bash
   # Increase Java heap size
   docker run -e JAVA_OPTS="-Xmx16g" ...
   ```

2. **Permission Issues**
   ```bash
   # Ensure output directory is writable
   chmod 755 output/
   ```

3. **VEP Cache Issues**
   ```bash
   # Verify VEP cache path is correctly mounted
   docker run --rm -v /path/to/vep:/vep:ro oligogenicity:latest ls /vep
   ```

### Getting Help

```bash
# Show general help
docker run --rm oligogenicity:latest --help

# Show command-specific help
docker run --rm oligogenicity:latest readvcfs --help
```

## Security Considerations

- The container runs as a non-root user (`appuser`) for security
- Mount data directories as read-only when possible
- Ensure sensitive data is properly secured on the host system
- Regular updates of the base image and dependencies 