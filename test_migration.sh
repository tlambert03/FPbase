#!/bin/bash
set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SOURCE_DB="fpbase_pre"
TARGET_DB="fpbase_migrated"

echo -e "${YELLOW}=== FPbase Migration Testing Script ===${NC}"
echo ""
echo "This script will:"
echo "  1. Drop '${TARGET_DB}' database if it exists"
echo "  2. Create a fresh copy from '${SOURCE_DB}'"
echo "  3. Run migration 0059 on the new database"
echo "  4. Leave '${SOURCE_DB}' completely untouched"
echo ""

# Check if source database exists
if ! psql -lqt | cut -d \| -f 1 | grep -qw "$SOURCE_DB"; then
    echo "ERROR: Source database '$SOURCE_DB' does not exist!"
    exit 1
fi

echo -e "${GREEN}Step 1: Dropping existing '${TARGET_DB}' database (if it exists)${NC}"
psql postgres -c "DROP DATABASE IF EXISTS ${TARGET_DB};" 2>/dev/null || true

echo -e "${GREEN}Step 2: Creating '${TARGET_DB}' as a copy of '${SOURCE_DB}'${NC}"
psql postgres -c "CREATE DATABASE ${TARGET_DB} WITH TEMPLATE ${SOURCE_DB};"

echo -e "${GREEN}Step 3: Running migration 0059 on '${TARGET_DB}'${NC}"

# Set DATABASE_URL to point to the new database
# This assumes PostgreSQL is running locally on default port
export DATABASE_URL="postgresql://localhost/${TARGET_DB}"

# Run the specific migration
uv run backend/manage.py migrate proteins 0059

echo ""
echo -e "${GREEN}=== Migration complete! ===${NC}"
echo ""
echo "Databases:"
echo "  - ${SOURCE_DB}: untouched (original test data)"
echo "  - ${TARGET_DB}: migrated (contains new schema)"
echo ""
echo "To inspect the migrated database:"
echo "  psql ${TARGET_DB}"
echo ""
echo "To run this script again:"
echo "  ./test_migration.sh"
