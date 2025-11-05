#!/bin/bash
set -e

echo "🚀 Setting up FPbase development environment..."

# Only run in remote Claude Code environment
if [ -z "$CLAUDE_CODE_REMOTE" ]; then
    exit 0
fi

# Fix sudo file ownership issues that can occur in cloud environments
echo "🔧 Fixing sudo permissions..."
chown root:root /etc/sudo.conf /etc/sudoers /etc/sudo_logsrvd.conf 2>/dev/null || true
chown -R root:root /etc/sudoers.d 2>/dev/null || true

# Install just command runner if not already installed
if ! command -v just &> /dev/null; then
    echo "📦 Installing just command runner..."
    cargo install just
fi

# Install PostgreSQL and development libraries (with network error tolerance)
echo "📦 Installing PostgreSQL and dependencies..."
apt-get update -qq 2>&1 | grep -v "^W:" || true
apt-get install -y --fix-missing \
    postgresql \
    postgresql-contrib \
    libpq-dev \
    gcc \
    python3-dev 2>&1 | grep -v "^W:" || true

# Install optional postgresql-server-dev-all (may fail due to dependencies)
apt-get install -y --fix-missing postgresql-server-dev-all 2>&1 | grep -v "^W:" || true

# Fix PostgreSQL SSL configuration (disable SSL for local testing)
echo "🔒 Configuring PostgreSQL SSL..."
if [ -f "/etc/postgresql/16/main/postgresql.conf" ]; then
    sed -i "s/^ssl = on/ssl = off/" /etc/postgresql/16/main/postgresql.conf
fi

# Fix PostgreSQL authentication (use trust for local connections)
echo "🔐 Configuring PostgreSQL authentication..."
if [ -f "/etc/postgresql/16/main/pg_hba.conf" ]; then
    # Change peer authentication to trust for local connections
    sed -i 's/^local\s\+all\s\+postgres\s\+peer/local   all             postgres                                trust/' /etc/postgresql/16/main/pg_hba.conf
    sed -i 's/^local\s\+all\s\+all\s\+peer/local   all             all                                     trust/' /etc/postgresql/16/main/pg_hba.conf
fi

# Fix file ownership for PostgreSQL
echo "📁 Fixing PostgreSQL file permissions..."
chown -R postgres:postgres /etc/postgresql/16/main/ 2>/dev/null || true
chown -R postgres:postgres /var/lib/postgresql/16/main/ 2>/dev/null || true
chown -R postgres:postgres /var/log/postgresql/ 2>/dev/null || true
chown postgres:postgres /var/run/postgresql/ 2>/dev/null || true

# Start PostgreSQL service
echo "🗄️  Starting PostgreSQL..."
service postgresql start

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
echo "🔧 Configuring PostgreSQL database..."
sudo -u postgres psql -c "ALTER USER postgres PASSWORD 'postgres';" 2>/dev/null || true
sudo -u postgres psql -c "CREATE DATABASE fpbase;" 2>/dev/null || echo "ℹ️  Database may already exist"

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
uv sync

# Install Node dependencies
echo "📦 Installing Node dependencies..."
pnpm install

# Install Playwright browsers for e2e tests (optional, can be slow)
echo "🎭 Installing Playwright browsers (this may take a while)..."
uv run playwright install chromium webkit 2>&1 | tail -5 || echo "⚠️  Warning: Playwright browser installation had issues, e2e tests may not work"

echo ""
echo "✅ Environment setup complete!"
echo ""
echo "You can now run:"
echo "  just test-py   # Run Python tests (165 tests)"
echo "  just test-js   # Run JavaScript tests (3 tests)"
echo "  just test-e2e  # Run end-to-end tests (49 tests, requires Playwright)"
echo "  just test      # Run all tests"
echo ""
