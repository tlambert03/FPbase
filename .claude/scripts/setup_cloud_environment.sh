#!/bin/bash
set -e

echo "🚀 Setting up FPbase development environment..."

# Only run in remote Claude Code environment
if [ -z "$CLAUDE_CODE_REMOTE" ]; then
    exit 0
fi

# Install PostgreSQL and development libraries
echo "📦 Installing PostgreSQL and dependencies..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
    postgresql \
    postgresql-contrib \
    postgresql-server-dev-all \
    libpq-dev \
    gcc \
    python3-dev

# Start PostgreSQL service
echo "🗄️  Starting PostgreSQL..."
sudo service postgresql start

# Wait for PostgreSQL to be ready
echo "⏳ Waiting for PostgreSQL to be ready..."
for i in {1..30}; do
    if sudo -u postgres psql -c '\q' 2>/dev/null; then
        echo "✅ PostgreSQL is ready"
        break
    fi
    if [ $i -eq 30 ]; then
        echo "❌ PostgreSQL failed to start"
        exit 1
    fi
    sleep 1
done

# Set password for postgres user and create database
echo "🔧 Configuring PostgreSQL..."
sudo -u postgres psql -c "ALTER USER postgres PASSWORD 'postgres';"
sudo -u postgres psql -c "CREATE DATABASE fpbase;" || echo "Database may already exist"

# Set environment variables
echo "🔑 Setting environment variables..."
if [ -n "$CLAUDE_ENV_FILE" ]; then
    cat >> "$CLAUDE_ENV_FILE" << 'EOF'
DATABASE_URL=postgres://postgres:postgres@localhost:5432/fpbase
UV_FROZEN=1
DJANGO_SETTINGS_MODULE=config.settings.test
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=fpbase
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
EOF
fi

# Install Python dependencies
echo "📚 Installing Python dependencies..."
uv sync --group test --no-dev

# Install Node dependencies
echo "📦 Installing Node dependencies..."
pnpm install

echo "✅ Environment setup complete!"
echo ""
echo "You can now run:"
echo "  uv run pytest                    # Run Python tests"
echo "  pnpm --filter @fpbase/spectra test:ci  # Run frontend tests"
